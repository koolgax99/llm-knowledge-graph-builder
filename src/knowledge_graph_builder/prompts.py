"""Centralized repository for all LLM prompts used in the knowledge graph system."""

# Phase 1: Main extraction prompts
MAIN_SYSTEM_PROMPT = """
You are an advanced AI system specialized in knowledge extraction and knowledge graph generation.
Your expertise includes identifying consistent entity references and meaningful relationships in text.
CRITICAL INSTRUCTION: All relationships (predicates) MUST be no more than 3 words maximum. Ideally 1-2 words. This is a hard limit.
"""

MAIN_USER_PROMPT = """
Your task: Read the text below (delimited by triple backticks) and identify all Subject-Predicate-Object (S-P-O) relationships in each sentence. Then produce a single JSON array of objects, each representing one triple.

Follow these rules carefully:

- Entity Consistency: Use consistent names for entities throughout the document. For example, if "Systemic Lupus Erythematosus" is mentioned as "sle", "lupus", and "rheumatic lupus" in different places, use a single consistent form (preferably the most complete one) in all triples.
- Atomic Terms: Identify distinct key terms (e.g., objects, diseases, locations, acronyms, people, conditions, concepts, feelings, environmental concepts). Avoid merging multiple ideas into one term (they should be as "atomistic" as possible).
- Unified References: Replace any pronouns (e.g., "he," "she," "it," "they," etc.) with the actual referenced entity, if identifiable.
- Pairwise Relationships: If multiple terms co-occur in the same sentence (or a short paragraph that makes them contextually related), create one triple for each pair that has a meaningful relationship.
- CRITICAL INSTRUCTION: Predicates MUST be 1-3 words maximum. Never more than 3 words. Keep them extremely concise.
- Ensure that all possible relationships are identified in the text and are captured in an S-P-O relation.
- Standardize terminology: If the same concept appears with slight variations (e.g., "artificial intelligence" and "AI"), use the most common or canonical form consistently.
- Make all the text of S-P-O text lower-case, even Names of people and places.
- If a person is mentioned by name, DO NOT create a relation to their location, profession and what they are known for (invented, wrote, started, title, etc.) if known.


Important Considerations:
- Aim for precision in entity naming - use specific forms that distinguish between similar but different entities
- Maximize connectedness by using identical entity names for the same concepts throughout the document
- Consider the entire context when identifying entity references
- ALL PREDICATES MUST BE 3 WORDS OR FEWER - this is a hard requirement

Output Requirements:

- Do not include any text or commentary outside of the JSON.
- Return only the JSON array, with each triple as an object containing "subject", "predicate", and "object".
- Make sure the JSON is valid and properly formatted.

Example of the desired output structure:

[
  {
    "subject": "Term A",
    "predicate": "relates to",  // Notice: only 2 words
    "object": "Term B"
  },
  {
    "subject": "Term C",
    "predicate": "uses",  // Notice: only 1 word
    "object": "Term D"
  }
]

Important: Only output the JSON array (with the S-P-O objects) and nothing else

Text to analyze (between triple backticks):
"""

MAIN_LUPUS_SYSTEM_PROMPT = """
You are an advanced AI system specialized in knowledge extraction and knowledge graph generation for lupus /sle / Systemic Lupus Erythematosus disease related.
Your expertise includes identifying consistent entity references and meaningful relationships in text related to the topic above.
CRITICAL INSTRUCTION: All relationships (predicates) MUST be no more than 3 words maximum. Ideally 1-2 words. This is a hard limit.
"""

MAIN_SCLC_SYSTEM_PROMPT = """
You are an advanced AI system specialized in knowledge extraction and knowledge graph generation for SCLC / Small Cell Lung Cancer / Neuroendocrines etc. related.
Your expertise includes identifying consistent entity references and meaningful relationships in text related to the topic above.
CRITICAL INSTRUCTION: All relationships (predicates) MUST be no more than 3 words maximum. Ideally 1-2 words. This is a hard limit.
"""

MAIN_LUPUS_USER_PROMPT = """
Your task: Read the research paper on Lupus (delimited by triple backticks) and identify only the meaningful Subject-Predicate-Object (S-P-O) relationships. Create a JSON array of these triplets while excluding author names, citations, publication dates, journal names, and other metadata.
Follow these rules carefully:

Focus on Medical Knowledge: Extract triplets related to lupus mechanisms, symptoms, treatments, diagnosis, and related medical concepts only.
Focus on Lupus and Food: Extract triplets related to lupus and how various foods, diet, chemicals affect progression.
Exclude Metadata: Do not create triplets about authors, institutions, publication details, or citation information.
Entity Consistency: Use consistent names for medical entities throughout. For example, if "Systemic Lupus Erythematosus" is also mentioned as "SLE" or just "Lupus", use the most complete form consistently in all triples.
Atomic Terms: Identify distinct medical terms (diseases, symptoms, organs, treatments, mechanisms, etc.). Keep terms atomistic rather than merging multiple concepts.
Unified References: Replace pronouns with the actual referenced medical entity.
CRITICAL INSTRUCTION: Predicates MUST be 1-3 words maximum. Never more than 3 words. Keep them extremely concise.
Standardize terminology: Use the canonical form of medical concepts consistently (e.g., always use "autoimmune disease" not "autoimmunity disorder" if referring to the same concept).
Make all the text of S-P-O triplets lowercase, even names of conditions and medical terms.
Focus on disease mechanisms, patient experiences, symptoms, treatments, organs affected, and medical processes rather than research methodology or study design.

Important Considerations:
- Aim for precision in entity naming - use specific forms that distinguish between similar but different entities
- Maximize connectedness by using identical entity names for the same concepts throughout the document
- Consider the entire context when identifying entity references
- ALL PREDICATES MUST BE 3 WORDS OR FEWER - this is a hard requirement

Output Requirements:

- Do not include any text or commentary outside of the JSON.
- Return only the JSON array, with each triple as an object containing "subject", "predicate", and "object".
- Make sure the JSON is valid and properly formatted.

Example of the desired output structure:
[
  {
    "subject": "lupus",
    "predicate": "affects",
    "object": "immune system"
  },
  {
    "subject": "hydroxychloroquine",
    "predicate": "treats",
    "object": "lupus symptoms"
  }
]
Important: Only output the JSON array with the medically meaningful S-P-O objects and nothing else. Exclude any triples related to paper authors, study methodology, publication details, or other scholarly metadata.
"""

MAIN_SCLC_USER_PROMPT = """
Your task: Read the research paper on Lupus (delimited by triple backticks) and identify only the meaningful Subject-Predicate-Object (S-P-O) relationships. Create a JSON array of these triplets while excluding author names, citations, publication dates, journal names, and other metadata.
Follow these rules carefully:

Focus on Medical Related Knowledge: Extract triplets related to SCLC / Neuroendocrines / Small cell lung cancer mechanisms, symptoms, treatments, diagnosis, and related medical concepts only. You can also extract triplets related to something that is indirectly related to SCLC, with third order or fourth order relationships.
Focus on SCLC and Food: Extract triplets related to SCLC and how various foods, diet, chemicals affect progression.
Exclude Metadata: Do not create triplets about authors, institutions, publication details, or citation information.
Entity Consistency: Use consistent names for medical entities throughout. For example, if "Small Cell Lung Cancer" is also mentioned as "SCLC", use the most complete form consistently in all triples.
Atomic Terms: Identify distinct medical terms (diseases, symptoms, organs, treatments, mechanisms, etc.). Keep terms atomistic rather than merging multiple concepts.
Unified References: Replace pronouns with the actual referenced medical entity.
CRITICAL INSTRUCTION: Predicates MUST be 1-3 words maximum. Never more than 3 words. Keep them extremely concise.
Standardize terminology: Use the canonical form of medical concepts consistently (e.g., always use "autoimmune disease" not "autoimmunity disorder" if referring to the same concept).
Make all the text of S-P-O triplets lowercase, even names of conditions and medical terms.
Focus on disease mechanisms, patient experiences, symptoms, treatments, organs affected, and medical processes rather than research methodology or study design.

Important Considerations:
- Aim for precision in entity naming - use specific forms that distinguish between similar but different entities
- Maximize connectedness by using identical entity names for the same concepts throughout the document
- Consider the entire context when identifying entity references
- ALL PREDICATES MUST BE 3 WORDS OR FEWER - this is a hard requirement

Output Requirements:

- Do not include any text or commentary outside of the JSON.
- Return only the JSON array, with each triple as an object containing "subject", "predicate", and "object".
- Make sure the JSON is valid and properly formatted.

Example of the desired output structure:
[
  {
    "subject": "lupus",
    "predicate": "affects",
    "object": "immune system"
  },
  {
    "subject": "hydroxychloroquine",
    "predicate": "treats",
    "object": "lupus symptoms"
  }
]
Important: Only output the JSON array with the medically meaningful S-P-O objects and nothing else. Exclude any triples related to paper authors, study methodology, publication details, or other scholarly metadata.
"""

MAIN_UV_USER_PROMPT = """
Your task: Analyze the UV forecasting research paper (delimited by triple backticks) and extract only the meaningful Subject-Predicate-Object (S-P-O) relationships related to UV forecasting. Present these as a JSON array of triplets.

Follow these specific guidelines:

1. Focus exclusively on UV knowledge: Extract triplets about UV radiation, UV forecasting, UVA/UVB, UV index, and directly related concepts.

2. Exclude all metadata: Do not create triplets about researchers, institutions, publication details, citations, or methodological information.

3. Entity consistency: Use uniform terminology for UV-related concepts throughout. If "ultraviolet radiation" appears as "UV radiation" or "UV," use the most complete form consistently.

4. Atomic terminology: Identify distinct UV-related terms (measurement methods, environmental factors, prediction models, atmospheric components, etc.). Keep terms atomistic rather than combining concepts.

5. Reference standardization: Replace pronouns with their actual referents.

6. Predicate brevity: CRITICAL - predicates MUST be 1-3 words maximum. Never exceed 3 words. Keep them extremely concise.

7. Terminology standardization: Use canonical forms of UV-related concepts consistently (e.g., always use "UV index" not "ultraviolet index" for the same concept).

8. Use lowercase for all text in the S-P-O triplets, even for proper terms and technical concepts.

9. Focus on UV forecasting mechanisms, measurement techniques, prediction models, atmospheric interactions, and environmental factors rather than study methodology.

Output requirements:
- Return only a properly formatted JSON array with no additional text or commentary
- Each triplet must be an object with "subject", "predicate", and "object" keys
- Ensure valid, properly formatted JSON syntax

Example of desired output structure:
[
  {
    "subject": "uv index",
    "predicate": "measures",
    "object": "radiation intensity"
  },
  {
    "subject": "ozone layer",
    "predicate": "filters",
    "object": "uv radiation"
  }
]

Important: Only output the JSON array with meaningful UV forecasting S-P-O relationships. Exclude any triplets related to paper authors, study design, publication details, or other scholarly metadata.
"""

# Phase 2: Entity standardization prompts
ENTITY_RESOLUTION_SYSTEM_PROMPT = """
You are an expert in entity resolution and knowledge representation.
Your task is to standardize entity names from a knowledge graph to ensure consistency.
"""


def get_entity_resolution_user_prompt(entity_list):
    return f"""
Below is a list of entity names extracted from a knowledge graph. 
Some may refer to the same real-world entities but with different wording.

Please identify groups of entities that refer to the same concept, and provide a standardized name for each group.
Return your answer as a JSON object where the keys are the standardized names and the values are arrays of all variant names that should map to that standard name.
Only include entities that have multiple variants or need standardization.

Entity list:
{entity_list}

Format your response as valid JSON like this:
{{
  "standardized name 1": ["variant 1", "variant 2"],
  "standardized name 2": ["variant 3", "variant 4", "variant 5"]
}}
"""


# Phase 3: Community relationship inference prompts
RELATIONSHIP_INFERENCE_SYSTEM_PROMPT = """
You are an expert in knowledge representation and inference. 
Your task is to infer plausible relationships between disconnected entities in a knowledge graph.
"""


def get_relationship_inference_user_prompt(entities1, entities2, triples_text):
    return f"""
I have a knowledge graph with two disconnected communities of entities. 

Community 1 entities: {entities1}
Community 2 entities: {entities2}

Here are some existing relationships involving these entities:
{triples_text}

Please infer 2-3 plausible relationships between entities from Community 1 and entities from Community 2.
Return your answer as a JSON array of triples in the following format:

[
  {{
    "subject": "entity from community 1",
    "predicate": "inferred relationship",
    "object": "entity from community 2"
  }},
  ...
]

Only include highly plausible relationships with clear predicates.
IMPORTANT: The inferred relationships (predicates) MUST be no more than 3 words maximum. Preferably 1-2 words. Never more than 3.
For predicates, use short phrases that clearly describe the relationship.
IMPORTANT: Make sure the subject and object are different entities - avoid self-references.
"""


# Phase 4: Within-community relationship inference prompts
WITHIN_COMMUNITY_INFERENCE_SYSTEM_PROMPT = """
You are an expert in knowledge representation and inference. 
Your task is to infer plausible relationships between semantically related entities that are not yet connected in a knowledge graph.
"""


def get_within_community_inference_user_prompt(pairs_text, triples_text):
    return f"""
I have a knowledge graph with several entities that appear to be semantically related but are not directly connected.

Here are some pairs of entities that might be related:
{pairs_text}

Here are some existing relationships involving these entities:
{triples_text}

Please infer plausible relationships between these disconnected pairs.
Return your answer as a JSON array of triples in the following format:

[
  {{
    "subject": "entity1",
    "predicate": "inferred relationship",
    "object": "entity2"
  }},
  ...
]

Only include highly plausible relationships with clear predicates.
IMPORTANT: The inferred relationships (predicates) MUST be no more than 3 words maximum. Preferably 1-2 words. Never more than 3.
IMPORTANT: Make sure that the subject and object are different entities - avoid self-references.
"""
