#Official ROOT Docker image
FROM rootproject/root-ubuntu16

USER root
RUN  apt-get update &&  apt-get install -y wget vim

# Run the following commands as super user (root):

#Instal FastJet
WORKDIR /
COPY installFastJet.sh .
RUN ./installFastJet.sh
ENV PATH "/fastjet/install/bin:${PATH}"
COPY installPythia.sh .
RUN ./installPythia.sh
ENV PATH "/pythia/install/bin:${PATH}"

RUN   apt-get install -y python-pip 
COPY requirements.txt .
RUN  sudo -H pip install --upgrade pip && sudo -H pip install --trusted-host pypi.python.org -r requirements.txt



ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password --gecos "Default user" \
            --uid ${NB_UID} ${NB_USER}

WORKDIR ${HOME}
USER ${NB_USER}
RUN  mkdir .jupyter && echo "c.NotebookApp.token = ''" > ${HOME}/.jupyter/jupyter_notebook_config.py
RUN  mkdir -p ${HOME}/plotter/eventFiles   ${HOME}/temp
#COPY --chown=jovyan exerciseNb ${HOME}/exerciseNb
#COPY --chown=jovyan exerciseNbExec ${HOME}/exerciseNbExec
#COPY --chown=jovyan exercisePy ${HOME}/exercisePy
EXPOSE 8888

# When starting the container and no command is started, run bash
#CMD ["/bin/bash"]
COPY --chown=jovyan plot.ipynb Makefile *.cpp  *.cc *.h *.hh ${HOME}/plotter/
COPY --chown=jovyan eventFiles/eventGG2H.txt ${HOME}/plotter/eventFiles
RUN cd plotter && make all

RUN jupyter nbextension enable --py widgetsnbextension
CMD ["jupyter", "notebook",  "--ip", "0.0.0.0"]
