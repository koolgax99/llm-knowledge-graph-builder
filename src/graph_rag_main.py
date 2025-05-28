import os
from neo4j_graph_rag.data_utils import load_kg_to_neo4j, add_embeddings_to_kg
from neo4j_graph_rag.knowledge_graph_rag import KnowledgeGraphRAG

# Configuration variables
file_name = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data_output/life_factor_lupus/life_factor_lupus_revised_prompt_2.json"  # Update this with your file path
neo4j_uri = "neo4j+s://8433b9d5.databases.neo4j.io"  # Update with your Neo4j URI
neo4j_user = "neo4j"  # Update with your Neo4j username
neo4j_password = "qEs51Ry7yZArHAfCcnTnLf914hoO-lHXjRV3VqtGKkI"  # Update with your Neo4j password
ollama_base_url = "http://149.165.169.239:11434/"  # Update with your Ollama base URL

# Alternatively, load from environment variables
# neo4j_password = os.environ.get("NEO4J_PASSWORD")
# openai_api_key = os.environ.get("OPENAI_API_KEY")

# Load knowledge graph to Neo4j
print(f"Loading knowledge graph from {file_name} into Neo4j...")
# load_kg_to_neo4j(
#     file_name,
#     neo4j_uri,
#     neo4j_user,
#     neo4j_password
# )
print("Knowledge graph loaded successfully.")

# Optional: Add embeddings to knowledge graph
add_embeddings = True  # Set to True if you want to add embeddings
if add_embeddings:
    print("Adding embeddings to knowledge graph nodes...")
    # add_embeddings_to_kg(
    #     neo4j_uri,
    #     neo4j_user,
    #     neo4j_password,
    #     ollama_base_url
    # )
    print("Embeddings added successfully.")
        
# Initialize the RAG system
print("Initializing RAG system...")
rag = KnowledgeGraphRAG(
    neo4j_uri,
    neo4j_user,
    neo4j_password,
    ollama_base_url,
    llm_model="gemma3:27b"
)
        
# Simple command-line interface for queries
print("\nKnowledge Graph RAG Query Interface")
print("Type 'exit' to quit")

while True:
    query = input("\nEnter your question: ")
    if query.lower() in ('exit', 'quit'):
        break
    
    try:
        answer = rag.answer_question(query)
        print(f"\nAnswer: {answer}")
    except Exception as e:
        print(f"\nError processing query: {e}")

