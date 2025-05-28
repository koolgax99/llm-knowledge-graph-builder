# Install required packages
# pip install neo4j langchain langchain_community langchain_ollama

import json
from neo4j import GraphDatabase
from langchain_ollama import OllamaEmbeddings

def load_kg_to_neo4j(json_file, uri, username, password):
    """Load knowledge graph from JSON file into Neo4j database."""
    
    # Load the knowledge graph data
    with open(json_file, 'r', encoding='utf-8') as f:
        triples = json.load(f)
    
    # Connect to Neo4j
    driver = GraphDatabase.driver(uri, auth=(username, password))
    
    # Create constraints for faster lookups (run once)
    with driver.session() as session:
        session.run("CREATE CONSTRAINT IF NOT EXISTS FOR (n:Entity) REQUIRE n.name IS UNIQUE")
    
    # Create nodes and relationships
    with driver.session() as session:
        for triple in triples:
            # Create subject and object nodes if they don't exist
            session.run("""
                MERGE (s:Entity {name: $subject})
                MERGE (o:Entity {name: $object})
                WITH s, o
                CREATE (s)-[r:RELATIONSHIP {type: $predicate}]->(o)
                """, 
                subject=triple["subject"], 
                object=triple["object"], 
                predicate=triple["predicate"]
            )
    
    driver.close()
    print(f"Loaded {len(triples)} triples into Neo4j")

def add_embeddings_to_kg(uri, username, password, ollama_base_url="http://localhost:11434"):
    """Add vector embeddings to entities and relationships using Ollama."""
    
    # Initialize Ollama embedding model
    # You can use models like 'nomic-embed-text' or other embedding models supported by Ollama
    embeddings = OllamaEmbeddings(
        model="nomic-embed-text",
        base_url=ollama_base_url
    )
    
    # Connect to Neo4j
    driver = GraphDatabase.driver(uri, auth=(username, password))
    
    # First, create a vector index (run once)
    with driver.session() as session:
        # Create vector index - adjust dimension based on your embedding model
        # nomic-embed-text uses 768 dimensions
        session.run("""
            CREATE VECTOR INDEX entity_embedding IF NOT EXISTS
            FOR (e:Entity) ON (e.embedding)
            OPTIONS {indexConfig: {
              `vector.dimensions`: 768,
              `vector.similarity_function`: 'cosine'
            }}
        """)
    
    # Get all entities and compute embeddings
    with driver.session() as session:
        result = session.run("MATCH (e:Entity) RETURN e.name AS name")
        entities = [record["name"] for record in result]
        
        # Process in batches 
        batch_size = 50
        for i in range(0, len(entities), batch_size):
            batch = entities[i:i+batch_size]
            
            # Compute embeddings
            entity_embeddings = embeddings.embed_documents(batch)
            
            # Update entities with embeddings
            for entity, embedding in zip(batch, entity_embeddings):
                session.run("""
                    MATCH (e:Entity {name: $name})
                    SET e.embedding = $embedding
                    """,
                    name=entity,
                    embedding=embedding
                )
            
            print(f"Processed {min(i+batch_size, len(entities))}/{len(entities)} entities")
    
    driver.close()