FROM python:3.11-slim  
  
# Set working directory  
WORKDIR /app  
  
# Copy requirements and install dependencies  
COPY requirements.txt .  
RUN pip install --no-cache-dir -r requirements.txt  
  
# Copy source code  
COPY src/ ./src/  
COPY pyproject.toml .  
COPY config.toml .   
  
# Create directories for data persistence  
RUN mkdir -p /app/data /app/output  

CMD ["bash"]