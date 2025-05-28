from langchain_ollama import OllamaEmbeddings, ChatOllama
from langchain_core.prompts import ChatPromptTemplate
from langchain_neo4j import Neo4jGraph, GraphCypherQAChain

class KnowledgeGraphRAG:
    def __init__(self, neo4j_uri, neo4j_username, neo4j_password, 
                 ollama_base_url="http://localhost:11434",
                 embedding_model="nomic-embed-text",
                 llm_model="llama3:latest"):
        """Initialize the Knowledge Graph RAG system."""
        
        # Set up Neo4j connection
        self.graph = Neo4jGraph(
            url=neo4j_uri,
            username=neo4j_username,
            password=neo4j_password
        )
        
        # Set up Ollama models
        self.embeddings = OllamaEmbeddings(
            model=embedding_model,
            base_url=ollama_base_url
        )
        
        self.llm = ChatOllama(
            model=llm_model,
            base_url=ollama_base_url,
            temperature=0
        )

        self.chain =  GraphCypherQAChain.from_llm(
            llm=self.llm,
            graph=self.graph,
            verbose=True,
            use_function_response=True,
            allow_dangerous_requests=True,
        )
    
    def vector_search(self, query, top_k=5):
        """Find entities related to the query using vector similarity."""
        
        # Embed the query
        query_embedding = self.embeddings.embed_query(query)
        
        # Search for similar entities
        result = self.graph.query("""
            CALL db.index.vector.queryNodes('entity_embedding', $top_k, $embedding)
            YIELD node, score
            RETURN node.name AS name, score
        """, {"top_k": top_k, "embedding": query_embedding})
        
        return [record for record in result]
    
    def expand_entities(self, entity_names, max_hops=2):
        """Explore relationships around the given entities."""
        
        # Construct parameters for query
        params = {"entities": entity_names, "max_hops": max_hops}
        
        # Get subgraph around entities
        result = self.graph.query("""
            MATCH path = (e:Entity)-[*1..$max_hops]-(related)
            WHERE e.name IN $entities
            RETURN path
        """, params)
        
        # Extract triples from the subgraph
        triples = []
        for record in result:
            path = record["path"]
            nodes = path.nodes
            rels = path.relationships
            
            for i in range(len(rels)):
                triple = {
                    "subject": nodes[i]["name"],
                    "predicate": rels[i].type,
                    "object": nodes[i+1]["name"]
                }
                if triple not in triples:
                    triples.append(triple)
        
        return triples
    
    def answer_question_old(self, question):
        """Generate an answer to the question based on knowledge graph information."""
        
        # Step 1: Find relevant entities using vector search
        relevant_entities = self.vector_search(question)
        entity_names = [e["name"] for e in relevant_entities]
        
        # Step 2: Expand to get related information
        knowledge_triples = self.expand_entities(entity_names)
        
        # Step 3: Format knowledge as context
        context = ""
        for triple in knowledge_triples:
            context += f"{triple['subject']} {triple['predicate']} {triple['object']}.\n"
        
        # Step 4: Generate answer using LLM
        prompt = ChatPromptTemplate.from_template("""
        You are an assistant that answers questions based on the provided knowledge.
        
        Knowledge:
        {context}
        
        Question: {question}
        
        Answer the question based only on the provided knowledge. If you cannot
        answer the question with the information provided, say "I don't have enough information to answer this question."
        """)
        
        chain = prompt | self.llm
        response = chain.invoke({"context": context, "question": question})
        
        return response.content
    
    def answer_question(self, question):
        response = self.chain.invoke({"query": question})
        print(response)
        return response['results']

