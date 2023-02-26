# Use Ubuntu 18.04 as the base image
FROM ubuntu:18.04

# Set the working directory to /app
WORKDIR /app

# Install Python 3 and pip
RUN apt-get update && \
    apt-get install -y python3 python3-pip bash

# Install the required packages
RUN pip3 install pytest mock pysam pyensembl

# Copy your code into the container
COPY your_code.py .

# Copy your test files into the container
COPY test_your_code.py .

# Echo the python version
CMD [ "python3", "-v" ]
