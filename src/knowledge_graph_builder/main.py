"""
Knowledge Graph Generator and Visualizer main module.
"""

import argparse
import json
import os
import sys

# # Add the parent directory to the Python path for imports
# sys.path.insert(
#     0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# )

from src.knowledge_graph_builder.config import load_config
from src.knowledge_graph_builder.llm import (
    call_llm,
    extract_json_from_text,
    call_ollama,
    call_openai,
)
from src.knowledge_graph_builder.visualization import (
    visualize_knowledge_graph,
    sample_data_visualization,
)
from src.knowledge_graph_builder.text_utils import chunk_text
from src.knowledge_graph_builder.entity_standardization import (
    standardize_entities,
    infer_relationships,
    limit_predicate_length,
)
from src.knowledge_graph_builder.prompts import (
    MAIN_SYSTEM_PROMPT,
    MAIN_USER_PROMPT,
    MAIN_LUPUS_USER_PROMPT,
    MAIN_UV_USER_PROMPT,
    MAIN_SCLC_USER_PROMPT,
    MAIN_SCLC_SYSTEM_PROMPT,
)


def process_with_llm(config, input_text, debug=False):
    """
    Process input text with LLM to extract triples.

    Args:
        config: Configuration dictionary
        input_text: Text to analyze
        debug: If True, print detailed debug information

    Returns:
        List of extracted triples or None if processing failed
    """
    # Use prompts from the prompts module
    system_prompt = MAIN_SCLC_SYSTEM_PROMPT
    user_prompt = MAIN_SCLC_USER_PROMPT
    user_prompt += f"```\n{input_text}```\n"

    # LLM configuration
    model = config["llm"]["model"]
    api_key = config["llm"]["api_key"]
    max_tokens = config["llm"]["max_tokens"]
    temperature = config["llm"]["temperature"]
    base_url = config["llm"]["base_url"]

    # Process with LLM
    metadata = {}
    response = call_openai(
        model, user_prompt, system_prompt, max_tokens, temperature, base_url
    )

    # Print raw response only if debug mode is on
    if debug:
        print("Raw LLM response:")
        print(response)
        print("\n---\n")

    # Extract JSON from the response
    result = extract_json_from_text(response)

    if result:
        # Validate and filter triples to ensure they have all required fields
        valid_triples = []
        invalid_count = 0

        for item in result:
            if (
                isinstance(item, dict)
                and "subject" in item
                and "predicate" in item
                and "object" in item
            ):
                # Add metadata to valid items
                valid_triples.append(dict(item, **metadata))
            else:
                invalid_count += 1

        if invalid_count > 0:
            print(
                f"Warning: Filtered out {invalid_count} invalid triples missing required fields"
            )

        if not valid_triples:
            print("Error: No valid triples found in LLM response")
            return None

        # Apply predicate length limit to all valid triples
        for triple in valid_triples:
            triple["predicate"] = limit_predicate_length(triple["predicate"])

        # Print extracted JSON only if debug mode is on
        if debug:
            print("Extracted JSON:")
            print(json.dumps(valid_triples, indent=2))  # Pretty print the JSON

        return valid_triples
    else:
        # Always print error messages even if debug is off
        print(
            "\n\nERROR ### Could not extract valid JSON from response: ",
            response,
            "\n\n",
        )
        return None


def process_text_in_chunks(config, full_text, debug=False):
    """
    Process a large text by breaking it into chunks with overlap,
    and then processing each chunk separately.

    Args:
        config: Configuration dictionary
        full_text: The complete text to process
        debug: If True, print detailed debug information

    Returns:
        List of all extracted triples from all chunks
    """
    # Get chunking parameters from config
    chunk_size = config.get("chunking", {}).get("chunk_size", 500)
    overlap = config.get("chunking", {}).get("overlap", 50)

    # Split text into chunks
    text_chunks = chunk_text(full_text, chunk_size, overlap)

    print("=" * 50)
    print("PHASE 1: INITIAL TRIPLE EXTRACTION")
    print("=" * 50)
    print(
        f"Processing text in {len(text_chunks)} chunks (size: {chunk_size} words, overlap: {overlap} words)"
    )

    # Process each chunk
    all_results = []
    for i, chunk in enumerate(text_chunks):
        print(
            f"Processing chunk {i + 1}/{len(text_chunks)} ({len(chunk.split())} words)"
        )

        # Process the chunk with LLM
        chunk_results = process_with_llm(config, chunk, debug)

        if chunk_results:
            # Add chunk information to each triple
            for item in chunk_results:
                item["chunk"] = i + 1

            # Add to overall results
            all_results.extend(chunk_results)
        else:
            print(f"Warning: Failed to extract triples from chunk {i + 1}")

    print(f"\nExtracted a total of {len(all_results)} triples from all chunks")

    # Apply entity standardization if enabled
    if config.get("standardization", {}).get("enabled", False):
        print("\n" + "=" * 50)
        print("PHASE 2: ENTITY STANDARDIZATION")
        print("=" * 50)
        print(
            f"Starting with {len(all_results)} triples and {len(get_unique_entities(all_results))} unique entities"
        )

        all_results = standardize_entities(all_results, config)

        print(
            f"After standardization: {len(all_results)} triples and {len(get_unique_entities(all_results))} unique entities"
        )

    # Apply relationship inference if enabled
    if config.get("inference", {}).get("enabled", False):
        print("\n" + "=" * 50)
        print("PHASE 3: RELATIONSHIP INFERENCE")
        print("=" * 50)
        print(f"Starting with {len(all_results)} triples")

        # Count existing relationships
        relationship_counts = {}
        for triple in all_results:
            relationship_counts[triple["predicate"]] = (
                relationship_counts.get(triple["predicate"], 0) + 1
            )

        print("Top 5 relationship types before inference:")
        for pred, count in sorted(
            relationship_counts.items(), key=lambda x: x[1], reverse=True
        )[:5]:
            print(f"  - {pred}: {count} occurrences")

        all_results = infer_relationships(all_results, config)

        # Count relationships after inference
        relationship_counts_after = {}
        for triple in all_results:
            relationship_counts_after[triple["predicate"]] = (
                relationship_counts_after.get(triple["predicate"], 0) + 1
            )

        print("\nTop 5 relationship types after inference:")
        for pred, count in sorted(
            relationship_counts_after.items(), key=lambda x: x[1], reverse=True
        )[:5]:
            print(f"  - {pred}: {count} occurrences")

        # Count inferred relationships
        inferred_count = sum(
            1 for triple in all_results if triple.get("inferred", False)
        )
        print(f"\nAdded {inferred_count} inferred relationships")
        print(f"Final knowledge graph: {len(all_results)} triples")

    return all_results


def get_unique_entities(triples):
    """
    Get the set of unique entities from the triples.

    Args:
        triples: List of triple dictionaries

    Returns:
        Set of unique entity names
    """
    entities = set()
    for triple in triples:
        if not isinstance(triple, dict):
            continue
        if "subject" in triple:
            entities.add(triple["subject"])
        if "object" in triple:
            entities.add(triple["object"])
    return entities


def process_file(config, input_file, output_file, debug=False):
    """
    Process a single file to extract triples and generate a knowledge graph.

    Args:
        config: Configuration dictionary
        input_file: Path to the input text file
        output_file: Path to the output JSON file
        debug: If True, print detailed debug information
    """
    try:
        with open(input_file, "r", encoding="utf-8") as f:
            input_text = f.read()
        
        print("=" * 50)
        print(f"Processing file: {input_file}")
        print("=" * 50)
    except Exception as e:
        print(f"Error reading input file {input_file}: {e}")
        return None

    # Process text in chunks
    result = process_text_in_chunks(config, input_text, debug)

    if result:
        # Save the raw data as JSON for potential reuse
        try:
            with open(output_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)

            # also save the raw data to txt file
            raw_output_file = output_file.replace(".json", ".txt")
            with open(raw_output_file, "w", encoding="utf-8") as f:
                for item in result:
                    # save result directly but without indentation for raw text
                    f.write(result + "\n")
            print(f"Saved processed triples to {output_file}")
            print(f"Saved raw knowledge graph data to {raw_output_file}")
        except Exception as e:
            print(f"Warning: Could not save raw data to {output_file}: {e}")
        return result
    else:
        print(f"Knowledge graph generation failed for file: {input_file}")
        return None


def main_kg_builder(
    input_path,
    output_dir,
    config_path="config.toml",
    experiment_name=None,
    debug=False,
    no_standardize=False,
    no_inference=False,
):
    """Main entry point for the knowledge graph generator.

    Args:
        input_path: Path to input text file or folder (required unless test is True)
        output_path: Output HTML file path
        config_path: Path to configuration file
        test: Generate a test visualization with sample data
        debug: Enable debug output (raw LLM responses and extracted JSON)
        no_standardize: Disable entity standardization
        no_inference: Disable relationship inference
    """
    # Load configuration
    config = load_config(config_path)
    if not config:
        print(f"Failed to load configuration from {config_path}. Exiting.")
        return

    # For normal processing, input path is required
    if not input_path:
        print("Error: input_path is required unless test is True")
        return

    # Override configuration settings with function arguments
    if no_standardize:
        config.setdefault("standardization", {})["enabled"] = False
    if no_inference:
        config.setdefault("inference", {})["enabled"] = False

    # Folder processing
    input_folder = input_path
    output_folder = os.path.join(output_dir, f"{experiment_name}_json_outputs")
    os.makedirs(output_folder, exist_ok=True)
    input_files = os.listdir(input_folder)

    for index, file_name in enumerate(input_files):
        
        print("=" * 50)
        print(f"Processing file {index + 1}/{len(input_files)}: {file_name}")
        print("=" * 50)

        if file_name.endswith(".txt"):
            input_file = os.path.join(input_folder, file_name)
            output_file = os.path.join(
                output_folder, f"{os.path.splitext(file_name)[0]}.json"
            )
            if output_file in os.listdir(output_folder):
                print(f"Skipping {output_file} as it already exists")
                continue
            
            result = process_file(config, input_file, output_file, debug)
            if result:
                # Visualize the combined knowledge graph
                output_file = output_file.replace(".json", ".html")
                stats = visualize_knowledge_graph(result, output_file, config=config)
                print("\nKnowledge Graph Statistics:")
                print(f"Nodes: {stats['nodes']}")
                print(f"Edges: {stats['edges']}")
                print(f"Communities: {stats['communities']}")
                print(
                    "\nTo view the visualization, open the following file in your browser:"
                )
                print(f"file://{os.path.abspath(output_file)}")
        else:
            print("No valid knowledge graphs were generated for this file")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Knowledge Graph Generator")
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to input folder containing text files",
    )
    parser.add_argument("--output", "-o", required=True, help="Output directory path")
    parser.add_argument(
        "--config", "-c", default="config.toml", help="Path to configuration file"
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", help="Enable debug output"
    )
    parser.add_argument(
        "--no-standardize", action="store_true", help="Disable entity standardization"
    )
    parser.add_argument(
        "--no-inference", action="store_true", help="Disable relationship inference"
    )
    args = parser.parse_args()

    main_kg_builder(
        args.input,
        args.output,
        config_path=args.config,
        debug=args.debug,
        no_standardize=args.no_standardize,
        no_inference=args.no_inference,
    )
