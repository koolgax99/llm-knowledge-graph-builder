import os
from PyPDF2 import PdfReader

def pdf_to_text(pdf_path, output_txt_path):
    """
    Converts a PDF file to a text file.

    :param pdf_path: Path to the input PDF file.
    :param output_txt_path: Path to the output text file.
    """
    if not os.path.exists(pdf_path):
        raise FileNotFoundError(f"The file {pdf_path} does not exist.")

    reader = PdfReader(pdf_path)
    with open(output_txt_path, 'w', encoding='utf-8') as txt_file:
        for page in reader.pages:
            txt_file.write(page.extract_text())
            txt_file.write("\n")  # Add a newline between pages

def process_pdf_folder(input_folder, output_folder):
    """
    Processes all PDF files in a folder and converts them to text files.

    :param input_folder: Path to the folder containing PDF files.
    :param output_folder: Path to the folder where text files will be saved.
    """
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"The folder {input_folder} does not exist.")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith(".pdf"):
            input_pdf_path = os.path.join(input_folder, file_name)
            file_name = file_name.replace(" ", "_")  # Replace spaces with underscores
            output_txt_path = os.path.join(output_folder, f"{os.path.splitext(file_name)[0]}.txt")
            try:
                print(f"Processing {input_pdf_path}...")
                pdf_to_text(input_pdf_path, output_txt_path)
                print(f"Saved text to {output_txt_path}")
            except Exception as e:
                print(f"Failed to process {input_pdf_path}: {e}")

if __name__ == "__main__":
    # Example usage
    input_path = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data_input_exp_4"  # Replace with your input folder or file path
    output_path = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data_input_exp_4/txt"  # Replace with your desired output folder path

    if os.path.isfile(input_path):
        # Single file processing
        output_file = os.path.join(output_path, f"{os.path.splitext(os.path.basename(input_path))[0]}.txt")
        try:
            pdf_to_text(input_path, output_file)
            print(f"PDF converted to text and saved at {output_file}")
        except Exception as e:
            print(f"An error occurred: {e}")
    elif os.path.isdir(input_path):
        # Folder processing
        try:
            process_pdf_folder(input_path, output_path)
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        print(f"The path {input_path} is neither a file nor a folder.")