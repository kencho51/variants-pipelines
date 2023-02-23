# Base image
FROM alpine:latest 

# Install necessary packages
RUN apk update && apk add \
        bash \
        curl \
        git \
        python3 \
        py3-pip

RUN pip3 install \
    pysam \
    vcf

CMD ["python3", "version"]