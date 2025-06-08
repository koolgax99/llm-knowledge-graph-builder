# from Bio import Entrez
# Entrez.email = "sanda.n@northeastern.edu"
# handle = Entrez.esearch(db="pubmed", retmax=10, term="opuntia[ORGN] accD", idtype="acc")
# record = Entrez.read(handle)
# handle.close()

# print(record['QueryTranslation'])

from marker.converters.pdf import PdfConverter
from marker.models import create_model_dict
from marker.config.parser import ConfigParser
from marker.output import text_from_rendered

config = {
    "output_format": "markdown",
    "disable_image_extraction": "true"
}
config_parser = ConfigParser(config)
from marker.models import create_model_dict

print(create_model_dict())

# converter = PdfConverter(
#     config=config_parser.generate_config_dict(),
#     artifact_dict=create_model_dict(),
#     processor_list=config_parser.get_processors(),
#     renderer=config_parser.get_renderer(),
#     llm_service=config_parser.get_llm_service()
# )
# rendered = converter("/home/exouser/masters-thesis/ai-knowledge-graph-main/data/input/data_input/pdf/lupus_and_climate/google_scholar/SLE Flares Seasonal Variations.pdf")
# text, meta, images = text_from_rendered(rendered)

# print("Metadata:", meta)

# # Save the markdown output to a file "text.md"
# with open("text.md", "w") as md_file:
#     md_file.write(text)


