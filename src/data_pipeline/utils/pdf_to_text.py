import os
import logging
from marker.converters.pdf import PdfConverter
from marker.models import create_model_dict
from marker.output import text_from_rendered
from marker.config.parser import ConfigParser


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pdf_converter.log'),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)


class PDFMarkerConverter:
    def __init__(self):
        """Initialize the converter and load models once."""
        logger.info("Loading marker models...")

        config = {
            "output_format": "markdown",
            "disable_image_extraction": "true"
        }
        config_parser = ConfigParser(config)

        self.converter = PdfConverter(
            config=config_parser.generate_config_dict(),
            artifact_dict=create_model_dict()
        )
        logger.info("Models loaded successfully!")

    def pdf_to_markdown(self, pdf_path, output_md_path):
        """
        Converts a PDF file to a markdown file using marker.

        :param pdf_path: Path to the input PDF file.
        :param output_md_path: Path to the output markdown file.
        """
        if not os.path.exists(pdf_path):
            raise FileNotFoundError(f"The file {pdf_path} does not exist.")

        try:
            # Convert PDF to markdown using marker
            rendered = self.converter(pdf_path)
            full_text, out_meta, images = text_from_rendered(rendered)

            # Write the markdown content to file
            with open(output_md_path, 'w', encoding='utf-8') as md_file:
                md_file.write(full_text)
            
            return out_meta
            
        except Exception as e:
            raise Exception(f"Failed to convert PDF using marker: {e}")

    def process_pdf_folder(self, input_folder, output_folder):
        """
        Processes all PDF files in a folder and converts them to markdown files.

        :param input_folder: Path to the folder containing PDF files.
        :param output_folder: Path to the folder where markdown files will be saved.
        """
        if not os.path.exists(input_folder):
            raise FileNotFoundError(f"The folder {input_folder} does not exist.")

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        processed_count = 0
        failed_count = 0

        for file_name in os.listdir(input_folder):
            if file_name.lower().endswith(".pdf"):
                input_pdf_path = os.path.join(input_folder, file_name)
                # Replace spaces with underscores and change extension to .md
                output_md_path = os.path.join(
                    output_folder, f"{os.path.splitext(file_name)[0]}.md")
                logger.info(f"Output markdown path: {output_md_path}")

                try:
                    logger.info(f"Processing {input_pdf_path}...")
                    metadata = self.pdf_to_markdown(input_pdf_path, output_md_path)
                    logger.info(f"✓ Saved markdown to {output_md_path}")
                    
                    processed_count += 1
                    
                except Exception as e:
                    logger.error(f"✗ Failed to process {input_pdf_path}: {e}")
                    failed_count += 1

        logger.info("Processing complete!")
        logger.info(f"Successfully processed: {processed_count} files")
        logger.info(f"Failed: {failed_count} files")


def main_pdf_converter(input_path, output_path):  
    """Main function to handle single file or folder processing."""
    
    # Initialize converter (models are loaded once)
    converter = PDFMarkerConverter()

    if os.path.isfile(input_path):
        # Single file processing
        if not input_path.lower().endswith('.pdf'):
            logger.error("Input file must be a PDF")
            return
            
        output_file = os.path.join(
            output_path, f"{os.path.splitext(os.path.basename(input_path))[0]}.md")
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        try:
            metadata = converter.pdf_to_markdown(input_path, output_file)
            logger.info(f"✓ PDF converted to markdown and saved at {output_file}")
            if metadata:
                logger.info(f"Metadata: {metadata}")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            
    elif os.path.isdir(input_path):
        # Folder processing
        try:
            converter.process_pdf_folder(input_path, output_path)
        except Exception as e:
            logger.error(f"An error occurred: {e}")
    else:
        logger.error(f"The path {input_path} is neither a file nor a folder.")


if __name__ == "__main__":

    # Configuration - Update these paths
    input_path = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/input/data_input/pdf/lupus_and_climate/google_scholar"
    output_path = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/input/data_input/txt/lupus_climate_google_scholar_marker"  # Changed to md_files
    main(input_path, output_path)