# Inspired from https://github.com/fbreitwieser/pavian

# Working with "release" codebase ----------------------------------------------
release-build:
	docker build \
	--no-cache \
	-f docker/release/Dockerfile \
	-t lianos/sparrow:release docker/release

release-run:
	docker run --rm -it -d --name sparrow-release \
	-v /Users/lianoglou/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v /Users/lianoglou/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v /Users/lianoglou/.config/rstudio:/home/rstudio/.config/rstudio \
	-v /Users/lianoglou/.gitconfig:/home/rstudio/.gitconfig \
	-v /Users/lianoglou/.gitignore_global:/Users/lianoglou/.gitignore_global \
	-p 8787:8787 \
	-p 8080:3838 \
	--entrypoint /init \
	--env SHINYPROXY_USERNAME=$(USER) \
	-e DISABLE_AUTH=true \
	lianos/sparrow:release && \
	sleep 5 && open http://localhost:8787

# Working with "devel" codebase ------------------------------------------------
devel-build:
	docker build \
	--no-cache \
	-f docker/devel/Dockerfile \
	-t lianos/sparrow:devel docker/devel

devel-run:
	docker run --rm -it -d --name sparrow-devel \
	-v ${HOME}/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v ${HOME}/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v ${HOME}/.config/rstudio:/home/rstudio/.config/rstudio \
	-v ${HOME}/.gitconfig:/home/rstudio/.gitconfig \
	-v ${HOME}/.gitignore_global:/Users/lianoglou/.gitignore_global \
	-p 8888:8787 \
	-p 8080:3838 \
	--entrypoint /init \
	--env SHINYPROXY_USERNAME=$(USER) \
	-e DISABLE_AUTH=true \
	lianos/sparrow:devel && \
	sleep 5 && open http://localhost:8787

devel-inspect:
	docker run --rm \
	-v ${HOME}/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v ${HOME}/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v ${HOME}/.config/rstudio:/home/rstudio/.config/rstudio \
	-v ${HOME}/.gitconfig:/home/rstudio/.gitconfig \
	-v ${HOME}/.gitignore_global:/Users/lianoglou/.gitignore_global \
	-i -t --entrypoint /bin/bash \
	lianos/sparrow:devel

