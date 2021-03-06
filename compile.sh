#!/bin/bash

echo "Set the GOPATH to '$GOPATH'"
export GOPATH=$(pwd)/Go
echo "Set the GOBIN to '$GOBIN'"
export GOBIN=$(pwd)/bin

if [[ ! -d "bin" ]]
then
    mkdir bin
else
    echo "./bin directory already exists"
fi

echo "Get Libaries"
#go get "github.com/go-gl/gl/v3.3-core/gl"
#go get "github.com/go-gl/glfw/v3.2/glfw"
#go get "github.com/go-gl/mathgl/mgl32"
#go get "github.com/fogleman/ln/ln"
go get "github.com/llgcode/draw2d/draw2dimg"
go get "github.com/paulmach/go.geo"

echo "Build Task"
go install Voronoi



