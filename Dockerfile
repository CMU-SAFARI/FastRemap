FROM ubuntu:21.10

RUN apt-get update -y --fix-missing
RUN apt-get upgrade -y

RUN apt-get install git make g++ zlib1g-dev libbz2-dev liblzma-dev -y

RUN mkdir /fastremap
RUN mkdir /input
RUN mkdir /output
WORKDIR /fastremap

RUN git clone --recursive https://github.com/CMU-SAFARI/FastRemap.git /fastremap
RUN cd zlib && ./configure && make && cd ..
RUN make
RUN cp FastRemap /usr/local/bin

RUN apt-get remove git -y
RUN apt-get autoremove -y
VOLUME /input
VOLUME /output
ENTRYPOINT ["/usr/local/bin/FastRemap"]

