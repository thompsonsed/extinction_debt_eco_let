# Note that there must be a tag
FROM thompsonsed/pycoalescence-circleci-0.0.1@sha256:a99f010094fecc71c5001ae9f9cdc79417bf3537ea49d8d98453f1d43dde4099
# install pip
RUN pip install --no-cache-dir notebook==5.*

# Create the user
ENV NB_USER extinction_debt_eco_let
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN pip install --no-cache-dir --user pycoalescence pandas

# Specify the default command to run
CMD jupyter notebook --ip 0.0.0.0 --notebook-dir=${HOME}
#CMD "jupyter", "notebook", "--ip", "0.0.0.0", "--notebook-dir=${HOME}"]