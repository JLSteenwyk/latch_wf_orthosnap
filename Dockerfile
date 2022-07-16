FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/wf-base:fbe8-main

# 
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

#create local directory to store output of cmd later
RUN mkdir /root/orthosnap_output/

COPY wf /root/wf

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN  sed -i 's/latch/wf/g' flytekit.config
RUN python3 -m pip install --upgrade latch
WORKDIR /root
