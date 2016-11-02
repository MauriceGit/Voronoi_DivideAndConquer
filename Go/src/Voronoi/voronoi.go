package main

import (
    . "Voronoi/HalfEdge"
    "fmt"
    "math/rand"
    "sort"
    //"math"
)

// Entry vertex for the whole voronoi diagram.
var voronoiStartVertex HEVertex

// I think, here we can dump our data and just keep references between them.
var heVertexList    []HEVertex
var heEdgeList      []HEEdge
var heFaceList      []HEFace

// There can only be 1 or 2 points in the list!!! This HAS to be guaranteed!
func createInitialVoronoi(points PointList) HEVertex {

    if len(points) == 1 {
        // A voronoi diagram with just one vertex is just this vertex with an infinite face around.
        return HEVertex{
            Pos:        points[0],
            ELeaving:   EmptyEdge,
        }
    } else {

    }

    return HEVertex{}
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func mergeVoronoi(left, right HEVertex) HEVertex {
    return HEVertex{}
}

// Voronoi divide and conquer entry point
func divideAndConquer(points PointList) HEVertex {

    // Recursion break on two points (what about 1 vs 3?)
    if len(points) <= 2 {
        return createInitialVoronoi(points)
    }

    // Split points in half
    halfLength := len(points)/2
    left  := divideAndConquer(points[:halfLength])
    right := divideAndConquer(points[halfLength:])

    return mergeVoronoi(left, right)

}

func main() {
    var r = rand.New(rand.NewSource(0))
    var pointList PointList

    for i:=0; i < 10; i++ {
        pointList = append(pointList, Vector{r.Float32()*100., r.Float32()*100.})
    }

    sort.Sort(pointList)

    voronoiStartVertex = divideAndConquer(pointList)

    fmt.Println(voronoiStartVertex)
}




