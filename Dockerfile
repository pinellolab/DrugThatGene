############################################################
# Dockerfile to build DrugThatGene
############################################################

# Set the base image to anaconda python 2.7
FROM continuumio/miniconda2:4.6.14

# File Author / Maintainer
MAINTAINER Luca Pinello 

ENV SHELL bash
#install dependencies
RUN conda install pip pandas==0.20.1
RUN pip install flask==0.12.2 gunicorn Flask-Excel
EXPOSE 9999


COPY DrugThatGene /DrugThatGene/DrugThatGene
COPY start_server_docker.sh /DrugThatGene/
COPY run.py /DrugThatGene/
WORKDIR DrugThatGene
CMD ["bash", "start_server_docker.sh"]
#CMD ["python","run.py"]
