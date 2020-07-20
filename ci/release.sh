#!/bin/bash
set -e
set -u

TAG=$(git describe)
CURRENT_VERSION=$(git describe | sed 's/^v//')
PREVIOUS_VERSION=$(git describe --abbrev=0 HEAD^)

# Build Docker image and tag
#
docker build -t guigolab/ggsashimi:$CURRENT_VERSION -f docker/Dockerfile .

# Run ci script in Docker container
#
docker run --entrypoint '/bin/bash' -v $PWD:$PWD -w $PWD guigolab/ggsashimi:$CURRENT_VERSION ci/run.sh

# Login to Docker Hub and push released versions
#
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
docker tag guigolab/ggsashimi:$CURRENT_VERSION guigolab/ggsashimi:latest
for t in latest $CURRENT_VERSION; do
    docker push "guigolab/ggsashimi:${t}"
done

# Get changelog and make release
#
git log --oneline --pretty="- %s" ...$PREVIOUS_VERSION | \
github-release release \
    --user guigolab \
    --repo ggsashimi \
    --tag  $TAG \
    --name "Version $CURRENT_VERSION" \
    --description - \
    --draft \
    --pre-release
