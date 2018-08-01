############################################################
# Dockerfile to build DrugThatGene
############################################################

# Set the base image to anaconda python 2.7
FROM continuumio/anaconda 

# File Author / Maintainer
MAINTAINER Luca Pinello 

ENV SHELL bash
#install dependencies
RUN pip install flask gunicorn Flask-Excel
EXPOSE 9999


COPY DrugThatGene /DrugThatGene/DrugThatGene
COPY start_server_docker.sh /DrugThatGene/
WORKDIR DrugThatGene
CMD ["bash", "start_server_docker.sh"]
