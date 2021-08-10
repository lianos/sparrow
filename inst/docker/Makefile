# Inspired from https://github.com/fbreitwieser/pavian

#release-build:
#	docker build -f release/Dockerfile \
#	--no-cache \
#	--build-arg GITHUB_PAT=$(GITHUB_PAT) \
#	-t lianos/sparrow:release .

release-build:
	docker build \
	--no-cache \
	-t lianos/sparrow:release release

devel-build:
	docker build \
	--no-cache \
	-t lianos/sparrow:devel devel

devel-run:
	docker run --rm -it -d --name sparrow-devel \
	-v /Users/lianoglou/workspace/Rpkgs/sparrow:/home/rstudio/sparrow \
	-v /Users/lianoglou/workspace/Rpkgs/sparrow.shiny:/home/rstudio/sparrow.shiny \
	-v /Users/lianoglou/.config/rstudio:/home/rstudio/.config/rstudio \
	-p 8787:8787 \
	-p 8080:3838 \
	--entrypoint /init \
	--env SHINYPROXY_USERNAME=$(USER) \
	-e DISABLE_AUTH=true \
	lianos/sparrow:devel && \
	sleep 5 && open http://localhost:8787

# devel-inspect:
#   docker run --rm \
# 	-v /Users/lianoglou/workspace/data/eyriedata:/datasets \
# 	-v /Users/lianoglou/workspace/facilebio/public/packages/FacileAnalysis:/home/rstudio/FacileAnalysis \
# 	-v /Users/lianoglou/workspace/facilebio/internal/DenaliEyrie:/home/rstudio/DenaliEyrie
# 	-i -t --entrypoint /bin/bash \
# 	-p 80:3838 $(TAG):$(VERSION)

