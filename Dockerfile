################## Base Image ##########
ARG PYTHON_VERSION="3.9.7"
FROM --platform=linux/amd64 python:${PYTHON_VERSION}

################## ARGUMENTS/Environments ##########
ARG BUILD_DATE
ARG BUILD_VERSION
ARG LICENSE="Apache-2.0"
ARG POSTPROCESSING_VARIANT_CALLS_VERSION="re-project-layout"
ARG VCS_REF

################## METADATA ########################
LABEL org.opencontainers.image.vendor="MSKCC"
LABEL org.opencontainers.image.authors="Eric Buehlere (buehlere@mskcc.org)"

LABEL org.opencontainers.image.created=${BUILD_DATE} \
    org.opencontainers.image.version=${BUILD_VERSION} \
    org.opencontainers.image.licenses=${LICENSE} \
    org.opencontainers.image.version.python=${PYTHON_VERSION} \ 
    org.opencontainers.image.vcs-ref=${VCS_REF} \ 
    org.opencontainers.image.version.conda=${CONDA_VERSION} \ 
    org.opencontainers.image.version.postprocessing_variant_calls=${POSTPROCESSING_VARIANT_CALLS_VERSION}

LABEL org.opencontainers.image.description="This container uses conda/conda/miniconda3 as the base image to build"

################## INSTALL ##########################
# download postprocessing_variant_calls 
RUN cd /opt \ 
    && git clone -b ${POSTPROCESSING_VARIANT_CALLS_VERSION} https://github.com/msk-access/postprocessing_variant_calls.git \ 
    && cd postprocessing_variant_calls 

# install postprocessing_variant_calls 
RUN cd /opt/postprocessing_variant_calls \
    && make deps-install

# by default /bin/sh
CMD ["/bin/sh"]