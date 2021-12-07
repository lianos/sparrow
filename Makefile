# Inspired from https://github.com/fbreitwieser/pavian

# Working with additional variables: https://stackoverflow.com/a/2826178/83761
# We can launch the devel version of the images like so:
# $ make run version=devel
version?=latest# `devel` or `latest`
rport?=8888# rstudio port
sport?=8080# shiny port

build:
	docker build \
	--no-cache \
	-f docker/${version}/Dockerfile \
	-t lianos/sparrow:${version} docker/${version}

run:
	docker run --rm -it -d --name sparrow-${version}-run \
	-v ${HOME}/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v ${HOME}/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v ${HOME}/.config/rstudio:/home/rstudio/.config/rstudio \
	-v ${HOME}/.gitconfig:/home/rstudio/.gitconfig \
	-v ${HOME}/.gitignore_global:/Users/lianoglou/.gitignore_global \
	-p ${rport}:8787 \
	-p ${sport}:3838 \
	--entrypoint /init \
	--env SHINYPROXY_USERNAME=$(USER) \
	-e DISABLE_AUTH=true \
	lianos/sparrow:${version} && \
	sleep 5 && open http://localhost:${rport}

inspect:
	docker run --rm --name sparrow-${version}-inspect \
	-v ${HOME}/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v ${HOME}/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v ${HOME}/.config/rstudio:/home/rstudio/.config/rstudio \
	-v ${HOME}/.gitconfig:/home/rstudio/.gitconfig \
	-v ${HOME}/.gitignore_global:/Users/lianoglou/.gitignore_global \
	-i -t --entrypoint /bin/bash \
	lianos/sparrow:${version}


