FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/wf-base:fbe8-main

# IQTree
RUN curl -L https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz -o iqtree-2.2.0-Linux.tar.gz &&\
    gunzip -cd iqtree-2.2.0-Linux.tar.gz | tar xfv - &&\
    mv iqtree-2.2.0-Linux/bin/iqtree2 /usr/bin &&\
    rm -rf iqtree-2.2.0-Linux iqtree-2.2.0-Linux.tar.gz

#create local directory to store output of cmd later
RUN mkdir /root/iqtree_output/

COPY wf /root/wf

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN  sed -i 's/latch/wf/g' flytekit.config
RUN python3 -m pip install --upgrade latch
WORKDIR /root
