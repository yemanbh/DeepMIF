FROM tensorflow/tensorflow:2.0.0-gpu-py3
MAINTAINER Name <YourEmail@abc.com>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# treated it like a Ubuntu OS:
#	update the apt-get;
#	set installation criteria of:
#		* y = yes to user prompts for proceeding (e.g. memory allocation);
#		* install only the recommended, not suggested, dependencies;
#		* include some essential installation packages for Ubuntu;
#		* python3-pip to allow pip be used to install the requirements
RUN apt-get update && apt-get install -y \
	--no-install-recommends \
	build-essential \
	python3-pip

# install packages for opencv
RUN apt-get update && apt install -y libsm6 \
                                      libxext6 \
                                      libxrender-dev
#################
# Install Python packages    #
#################
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

RUN mkdir code
COPY . code

ENTRYPOINT ["python", "code/run_deep_learning_pipeline.py"]

