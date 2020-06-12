FROM php:7-apache

MAINTAINER max@maxhebditch.co.uk

COPY ./html/ /var/www/html/
ADD ./code /var/www/html/code/
ADD ./html/css /var/www/html/css/

RUN apt-get -y update
RUN apt-get -y install python3 python3-pip zip
RUN cd /var/www/html/code && pip3 install -r requirements.txt

RUN mkdir /var/www/html/results
RUN chmod a+rwx -R /var/www/html/results
