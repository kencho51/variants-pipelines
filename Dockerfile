# Use Ubuntu 18.04 as the base image
FROM ubuntu:18.04

# Set the working directory to /app
#WORKDIR /app

# Install Python 3 and pip
RUN apt-get update && \ 
    apt-get install -y python3 python3-pip bash zlib1g-dev libbz2-dev liblzma-dev

# Upgrade pip3 and install polars
RUN pip3 install --upgrade pip 
RUN pip3 install polars

# Install pandas first
RUN pip3 install pandas==1.1.5

# Install the latest pyensembl
RUN pip3 install pyensembl==1.8.3


# Install the required packages
RUN pip3 install pytest mock pysam

# Copy your code into the container
#dockerCOPY . /app

# Echo the python version
CMD [ "python3", "-v" ]
