package main

import (
    . "Voronoi/Vector"
    . "Voronoi/HalfEdge"
    "fmt"
    "math"
    "os"
    "math/rand"
    "sort"
    //"math"
    "github.com/llgcode/draw2d/draw2dimg"
    "image"
    "image/color"
    "image/png"
    "time"
    "strconv"
    "errors"
    "container/list"
)

// HAS to be an infinite face!
type VoronoiEntryFace FaceIndex

// List of indides to the faces, that form the convex hull.
type ConvexHull []FaceIndex

type Voronoi struct {
    vertices           []HEVertex
    firstFreeVertexPos VertexIndex
    edges              []HEEdge
    firstFreeEdgePos   EdgeIndex
    faces              []HEFace
    firstFreeFacePos   FaceIndex
}

type ChainElem struct {
    intersection    Vector
    edgeP           EdgeIndex
    edgeQ           EdgeIndex
    p               FaceIndex
    q               FaceIndex
    bisector        Edge
}


var g_recursions int = 0
var g_drawImages bool = true;
var g_pprint bool = true;
var g_freeEdgePositions = list.New()
var g_freeVertexPositions = list.New()

//var emptyVector = Vector{-1, -1, -1}
//var emptyF = HEFace{emptyVector, -1}
//var emptyE = HEEdge{-1, -1, -1, -1, -1, Edge{emptyVector, emptyVector}}
//var emptyV = HEVertex{emptyVector}

var emptyF = HEFace{}
var emptyE = HEEdge{}
var emptyV = HEVertex{}

////////////////////////////////////////////////////////////////////////
//  Pretty Print the  Voronoi Attributes
////////////////////////////////////////////////////////////////////////

func (v *Voronoi) pprint () {
    if !g_pprint {
        return
    }
    fmt.Println("Voronoi:")
    fmt.Printf("   Vertices (%v):\n", v.firstFreeVertexPos)

    for i,ve := range v.vertices {
        if ve != emptyV {
            fmt.Printf("\t%v:\tPos: %v\n", i, ve.Pos)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Edges    (%v):\n", v.firstFreeEdgePos)
    for i,e := range v.edges {
        if e != emptyE {
            fmt.Printf("\t%v:\tOrigin: %v,\tTwin: %v,\tPrev: %v,\tNext: %v,\tFace: %v\n", i, e.VOrigin, e.ETwin, e.EPrev, e.ENext, e.FFace)
        }
    }
    fmt.Printf("\n")

    fmt.Printf("   Faces    (%v):\n", v.firstFreeFacePos)
    for i,f := range v.faces {
        if f != emptyF {
            fmt.Printf("\t%v:\tRefPoint: %v,\tEdge:%v\n", i, f.ReferencePoint, f.EEdge)
        }
    }
    fmt.Printf("\n")
}

////////////////////////////////////////////////////////////////////////
//  Draw an Image of the Voronoi Tessellation
////////////////////////////////////////////////////////////////////////

type Circle struct {
    X, Y, R float64
}

func (c *Circle) insideCircle(x, y float64) bool {
    var dx, dy float64 = c.X - x, c.Y - y
    d := math.Sqrt(dx*dx+dy*dy) / c.R
    return d <= 1
}

func drawCircle(m *image.RGBA, posX, posY, radius int, c color.RGBA) {
    cr := &Circle{float64(posX), float64(posY), float64(radius)}

    for x := posX-radius; x < posX+radius; x++ {
        for y := posY-radius; y < posY+radius; y++ {
            if cr.insideCircle(float64(x), float64(y)) {
                m.Set(x, y, c)
            }
        }
    }
}

func (v *Voronoi)createImage(filename string, whole bool) {

    if !g_drawImages {
        return
    }

    var w, h int = 1000, 1000
    scale := 10.0

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // Edges
    // yellow
    //c := color.RGBA{255,255,0,255}
    // blue
    c := color.RGBA{0,0,255,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(1)
    gc.SetStrokeColor(c)

    for i,e := range v.edges {
        if e != emptyE {
            e2   := v.edges[e.ETwin]
            if e2 != emptyE {

                edge := Edge{}
                v1   := e.VOrigin
                v2   := e2.VOrigin
                switch {
                    // Best case. We have both endpoints. And none of them Infinity ones.
                    case v1.Valid() && v2.Valid() && v.vertices[v1].Pos != InfinitePoint && v.vertices[v2].Pos != InfinitePoint:
                        //fmt.Printf("Edge in question: \n\t%v\n\t%v\n\t%v, %v\n", e, e2, v1, v2)
                        edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                    // We have the "left" endpoint.
                    case v1.Valid() && v.vertices[v1].Pos != InfinitePoint && !v2.Valid():
                        edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                    // We have the "right" endpoint.
                    case !v1.Valid() && v2.Valid() && v.vertices[v2].Pos != InfinitePoint:
                        edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                    // We don't have any endpoints.
                    default:
                        //i = i
                        // amplified line. Exceeding all boundaries. Infinite line.
                        edge = createLine(v, EdgeIndex(i), false)
                }



                if edge == (Edge{}) {
                    //continue
                }

                if i == 50 || i == 51 {
                    //fmt.Printf("coloring now for something\n")
                    //tmpC := color.RGBA{0,0,0,255}
                    //gc.SetStrokeColor(tmpC)
                } else {
                    //gc.SetStrokeColor(c)
                }

                //fmt.Printf("From %v|%v ----> %v|%v\n", edge.Pos.X*scale., float64(h) - edge.Pos.Y*scale., Add(edge.Pos, edge.Dir).X*scale., float64(h) - Add(edge.Pos, edge.Dir).Y*scale.)

                gc.MoveTo(edge.Pos.X*scale, float64(h) - edge.Pos.Y*scale)
                gc.LineTo(Add(edge.Pos, edge.Dir).X*scale, float64(h) - Add(edge.Pos, edge.Dir).Y*scale)
                gc.FillStroke()
                gc.Close()
            }
        }
    }

    // Faces/Reference Points!
    normalC := color.RGBA{255,0,0,255}
    //tmpC := color.RGBA{0,0,0,255}
    c = normalC
    for i,f := range v.faces {
        if i == 15 {
            //fmt.Printf("coloring face now for something\n")
            //c = tmpC
        } else {
            //c = normalC
        }
        drawCircle(m, int(f.ReferencePoint.X*scale), h-int(f.ReferencePoint.Y*scale), 2, c)
    }

    // Vertices between edges
    c = color.RGBA{0,255,0,255}
    for _,ve := range v.vertices {
        if ve != emptyV {
            drawCircle(m, int(ve.Pos.X*scale), h-int(ve.Pos.Y*scale), 2, c)
        }
    }

    f, err := os.OpenFile(filename + ".png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}


func drawEdgeList(v *Voronoi, edges []EdgeIndex, filename string) {
    if !g_drawImages {
        return
    }

    var w, h int = 1000, 1000

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // Edges
    // yellow
    //c := color.RGBA{255,255,0,255}
    // blue
    c := color.RGBA{0,100,100,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(2)
    gc.SetStrokeColor(c)

    for i,eI := range edges {

        if eI == EmptyEdge {
            continue
        }

        e := v.edges[eI]

        //var tmp HEEdge
        tmp := HEEdge{-1, -1, -1, -1, -1, Edge{Vector{-1, -1, -1}, Vector{-1, -1, -1}}}
        e2   := v.edges[e.ETwin]
        if e != tmp && e2 != tmp {

            edge := Edge{}
            v1   := e.VOrigin
            v2   := e2.VOrigin
            switch {
                // Best case. We have both endpoints. And none of them Infinity ones.
                case v1.Valid() && v2.Valid() && v.vertices[v1].Pos != InfinitePoint && v.vertices[v2].Pos != InfinitePoint:
                    //fmt.Printf("Edge in question: \n\t%v\n\t%v\n\t%v, %v\n", e, e2, v1, v2)
                    edge = Edge{v.vertices[v1].Pos, Sub(v.vertices[v1].Pos, v.vertices[v2].Pos)}

                // We have the "left" endpoint.
                case v1.Valid() && v.vertices[v1].Pos != InfinitePoint && !v2.Valid():
                    edge = Edge{v.vertices[v1].Pos, e.TmpEdge.Dir}

                // We have the "right" endpoint.
                case !v1.Valid() && v2.Valid() && v.vertices[v2].Pos != InfinitePoint:
                    edge = Edge{v.vertices[v2].Pos, e2.TmpEdge.Dir}

                // We don't have any endpoints.
                default:
                    //i = i
                    // amplified line. Exceeding all boundaries. Infinite line.
                    edge = createLine(v, EdgeIndex(i), false)
            }



            if edge == (Edge{}) {
                //continue
            }

            //fmt.Printf("From %v|%v ----> %v|%v\n", edge.Pos.X*10., float64(h) - edge.Pos.Y*10., Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)

            gc.MoveTo(edge.Pos.X*10., float64(h) - edge.Pos.Y*10.)
            gc.LineTo(Add(edge.Pos, edge.Dir).X*10., float64(h) - Add(edge.Pos, edge.Dir).Y*10.)
            gc.FillStroke()
            gc.Close()
        }
    }

    f, err := os.OpenFile(filename + ".png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

func (v *Voronoi)drawFaces() {
    if !g_drawImages {
        return
    }

    for i,f := range v.faces {
        var edgeIndexList []EdgeIndex

        startEdge := f.EEdge

        edgeIndexList = append(edgeIndexList, startEdge)

        edgeCount := 0
        e := v.edges[startEdge].ENext

        for e != EmptyEdge && e != startEdge {
            edgeIndexList = append(edgeIndexList, e)
            if edgeCount > 50 {
                fmt.Printf("error because of too large edge count?")
                break
            }

            e = v.edges[e].ENext
            edgeCount += 1
        }

        drawEdgeList(v, edgeIndexList, "edges_face_" + fmt.Sprintf("%v", i))

    }

}

func drawDividingChain(chain []ChainElem) {
    if !g_drawImages {
        return
    }

    var w, h int = 1000, 1000
    scaleX := 10.0
    scaleY := 10.0

    m := image.NewRGBA(image.Rect(0, 0, w, h))

    // blue
    c := color.RGBA{255,0,0,255}
    // green
    //c := color.RGBA{0,255,0,255}
    gc := draw2dimg.NewGraphicContext(m)

    gc.SetLineWidth(2)

    gc.SetStrokeColor(c)

    empty := Vector{}
    for i,ch := range chain {
        if ch.intersection != empty {
            drawCircle(m, int(ch.intersection.X*scaleX), h-int(ch.intersection.Y*scaleY), 3, c)

            if i > 0 {
                gc.MoveTo(chain[i-1].intersection.X*scaleX, float64(h) - chain[i-1].intersection.Y*scaleY)
                gc.LineTo(ch.intersection.X*scaleX, float64(h) - ch.intersection.Y*scaleY)
                gc.FillStroke()
                gc.Close()
            }
        }

    }

    f, err := os.OpenFile("dividing_chain.png", os.O_WRONLY|os.O_CREATE, 0600)
    if err != nil {
        fmt.Println(err)
        return
    }
    defer f.Close()
    png.Encode(f, m)
}

// Verifies that the Voronoi is not corrupted or has miscalculated edges/faces
// or any other unvalid stuff.
func (v *Voronoi) Verify() error {

    for i,e := range v.edges {
        if e != emptyE {

            // Every valid edge MUST have a valid Twin edge!
            if e.ETwin == EmptyEdge || v.edges[e.ETwin] == emptyE {
                return errors.New(fmt.Sprintf("Edge %v: %v has an invalid twin edge", i, e))
            }

            // Twins must refer to each other!
            if i != int(v.edges[e.ETwin].ETwin) {
                return errors.New(fmt.Sprintf("Edge %v and his Twin %v don't refer to each other", i, e.ETwin))
            }

            // Check if the origin vertex is valid (if there is one)
            if e.VOrigin != EmptyVertex && e.VOrigin != InfiniteVertex && v.vertices[e.VOrigin] == emptyV {
                return errors.New(fmt.Sprintf("Origin vertex of edge %v: %v is invalid", i, e))
            }

            // Check, if the face is valid
            if e.FFace == EmptyFace || v.faces[e.FFace] == emptyF {
                return errors.New(fmt.Sprintf("Face of edge %v: %v is invalid", i, e))
            }

            // if VOrigin is Empty, it MUST be the referenced edge for the face of the edge!
            if (e.VOrigin == EmptyVertex || e.VOrigin == InfiniteVertex) && i != int(v.faces[e.FFace].EEdge) {
                return errors.New(fmt.Sprintf("A first edge %v: %v must be referenced as first edge by the corresponding face %v!", e, i, e.FFace))
            }

            // Check, if the next edge is valid
            if e.ENext != EmptyEdge && v.edges[e.ENext] == emptyE {
                return errors.New(fmt.Sprintf("Next edge for edge %v: %v is invalid", i, e))
            }

            // Check, if the prev edge is valid
            if e.EPrev != EmptyEdge && v.edges[e.EPrev] == emptyE {
                return errors.New(fmt.Sprintf("Prev edge for edge %v: %v is invalid", i, e))
            }

            // If the prev edge is valid
            if e.VOrigin != EmptyVertex && e.VOrigin != InfiniteVertex && e.EPrev == EmptyEdge {
                return errors.New(fmt.Sprintf("Prev edge for edge %v: %v is not defined but should be", i, e))
            }

            // If this->next is not empty, the next one must have a previous edge.
            if e.ENext != EmptyEdge && v.edges[e.ENext].EPrev == EmptyEdge {
                return errors.New(fmt.Sprintf("Next edge for edge %v must have a non-empty prev", i))
            }

            // The this edge must correspond to the next->prev edge
            if e.ENext != EmptyEdge && v.edges[e.ENext].EPrev != EmptyEdge && v.edges[e.ENext].EPrev != EdgeIndex(i) {
                return errors.New(fmt.Sprintf("The edge %v has %v as next edge, which has %v as his previous. They must reference each other!", i, e.ENext, v.edges[e.ENext].EPrev))
            }



        }
    }

    for i,f := range v.faces {
        if f != emptyF {

            // Check for valid reference point
            if f.ReferencePoint == InfinitePoint {
                return errors.New(fmt.Sprintf("Face %v: %v has infinite reference point", i, f))
            }

            // Check for valid EEdge
            if f.EEdge == EmptyEdge || v.edges[f.EEdge] == emptyE {
                return errors.New(fmt.Sprintf("Face %v: %v points to invalid edge", i, f))
            }

            // If the edge is the first one (no/infinite start-vertex), it's all good.
            // Otherwise check, that the edges actually go all the way around the face!
            if v.edges[f.EEdge].VOrigin != EmptyVertex && v.edges[f.EEdge].VOrigin != InfiniteVertex {
                startEdge := f.EEdge
                edgeCount := 0
                e := v.edges[startEdge].ENext

                for e != EmptyEdge && e != startEdge {

                    if int(v.edges[e].FFace) != i {
                        return errors.New(fmt.Sprintf("Edge %v of the face %v: %v does not point to the right face!", e, i, f))
                    }

                    if edgeCount > 50 {
                        return errors.New(fmt.Sprintf("Looping around the edges of the face %v: %v most likely does not stop (infinite loop)", i, f))
                    }

                    e = v.edges[e].ENext
                    edgeCount += 1
                }
            }
        }
    }

    return nil
}

////////////////////////////////////////////////////////////////////////
//  Calculate a Voronoi Tessellation
////////////////////////////////////////////////////////////////////////

// Calculates a convex hull for the given voronoi. Runs in O(n).
// The convex hull consists of a list of Reference points.
// As there is an exactly 1:1 representation of Reference points and HEFaces,
// a list of face indices is returned in Order (counter-clockwise) of the convex hull.
func (v *Voronoi)ConvexHull(fEntry VoronoiEntryFace) ConvexHull {

    startFace := FaceIndex(fEntry)
    convexHullList := []FaceIndex{FaceIndex(startFace)}

    // If the voronoi consists of just one refPoint, there is no edge.
    if v.faces[startFace].EEdge != EmptyEdge {

        nextFace := v.edges[v.edges[v.faces[startFace].EEdge].ETwin].FFace

        convexHullList = append(convexHullList, nextFace)

        for {
            nextFace = v.edges[v.edges[v.faces[nextFace].EEdge].ETwin].FFace

            // Once all around.
            if nextFace == startFace {
                return convexHullList
            }

            lastFace := convexHullList[len(convexHullList)-2]
            // Special case of linear voronoi.
            if nextFace == lastFace {
                return convexHullList
            }

            convexHullList = append(convexHullList, nextFace)
        }
    }
    return convexHullList
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createFace(refPoint Vector, eEdge EdgeIndex) FaceIndex {
    v.faces[v.firstFreeFacePos] = HEFace {
        ReferencePoint: refPoint,
        EEdge:          eEdge,
    }
    v.firstFreeFacePos += 1
    return v.firstFreeFacePos-1
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createEdge(vOrigin VertexIndex, eTwin, ePrev, eNext EdgeIndex, fFace FaceIndex, tmpEdge Edge) EdgeIndex {

    index := v.firstFreeEdgePos
    if g_freeEdgePositions.Len() > 0 {
        elem := g_freeEdgePositions.Front()
        index = elem.Value.(EdgeIndex)
        // Now, this has to be a O(1) operations, because we Always pop the very first element of the list!
        // If the removal of the first element turns out to be O(n), we have to change to a better linked
        // list implementation!
        g_freeEdgePositions.Remove(elem)
    } else {
        v.firstFreeEdgePos += 1
    }

    v.edges[index] = HEEdge {
        VOrigin:    vOrigin,
        ETwin:      eTwin,
        ENext:      eNext,
        EPrev:      ePrev,
        FFace:      fFace,
        TmpEdge:    tmpEdge,
    }
    return index
}

func (v *Voronoi)deleteEdgePair(e1, e2 EdgeIndex) {
    fmt.Printf("DELETE E: %v | %v\n", e1, e2)

    // Delete all references to the edges.
    if v.edges[e1].EPrev != EmptyEdge && v.edges[v.edges[e1].EPrev].ENext == e1 {
       v.edges[v.edges[e1].EPrev].ENext = EmptyEdge
    }
    if v.edges[e1].ENext != EmptyEdge && v.edges[v.edges[e1].ENext].EPrev == e1 {
       v.edges[v.edges[e1].ENext].EPrev = EmptyEdge
    }
    if v.edges[e2].EPrev != EmptyEdge && v.edges[v.edges[e2].EPrev].ENext == e2 {
       v.edges[v.edges[e2].EPrev].ENext = EmptyEdge
    }
    if v.edges[e2].ENext != EmptyEdge && v.edges[v.edges[e2].ENext].EPrev == e2 {
       v.edges[v.edges[e2].ENext].EPrev = EmptyEdge
    }

    // Delete the edges themself.
    v.edges[e1]  = emptyE
    v.edges[e2]  = emptyE
    g_freeEdgePositions.PushBack(e1)
    g_freeEdgePositions.PushBack(e2)
}

func (v *Voronoi)deleteVertex(ve VertexIndex) {
    fmt.Printf("DELETE V: %v\n", ve)
    v.vertices[ve] = emptyV

    g_freeVertexPositions.PushBack(ve)
}

// So we can get a pointer of some data structure? What about scope issues?
func (v *Voronoi)createVertex(pos Vector) VertexIndex {
    index := v.firstFreeVertexPos
    if g_freeVertexPositions.Len() > 0 {
        elem := g_freeVertexPositions.Front()
        index = elem.Value.(VertexIndex)
        // Now, this has to be a O(1) operations, because we Always pop the very first element of the list!
        // If the removal of the first element turns out to be O(n), we have to change to a better linked
        // list implementation!
        g_freeVertexPositions.Remove(elem)
    } else {
        v.firstFreeVertexPos += 1
    }

    v.vertices[index] = HEVertex {
        Pos:        pos,
    }

    return index
}

// Calculates the face with the 'best' (y-coord) reference point.
// Can be used to calculate the lowest or highest point with the appropriate function.
// Returns the index of the convex hull list called upon.
func (ch ConvexHull) bestFaceIndex(v *Voronoi, isBetter func(v1, v2 Vector) bool) int {
    bestFace := 0
    for i,_ := range ch {
        if isBetter(v.faces[ch[i]].ReferencePoint, v.faces[ch[bestFace]].ReferencePoint) {
            bestFace = i
        }
    }
    return bestFace
}

// Creates a line from the HEEdge. Depeding on its state from the edge or the TmpEdge.
func createLine (v *Voronoi, e EdgeIndex, amplified bool) Edge {
    switch {
        case !v.edges[e].VOrigin.Valid() && v.edges[v.edges[e].ETwin].VOrigin.Valid():
            edge := v.edges[e].TmpEdge
            edge.Dir = Mult(edge.Dir, 50)
            return Amplify(edge, 50)
        case v.edges[e].VOrigin.Valid() && !v.edges[v.edges[e].ETwin].VOrigin.Valid():
            edge := v.edges[e].TmpEdge
            edge.Dir = Mult(edge.Dir, 50)
            return Amplify(edge, 50)
        case !v.edges[e].VOrigin.Valid() && !v.edges[v.edges[e].ETwin].VOrigin.Valid():
             return Amplify(v.edges[e].TmpEdge, 50.0)
        default:
            return Edge {
                Pos: v.vertices[v.edges[e].VOrigin].Pos,
                Dir: Sub(v.vertices[v.edges[v.edges[e].ETwin].VOrigin].Pos, v.vertices[v.edges[e].VOrigin].Pos),
            }
    }
}

// Calculates the highest intersection of the given bisector and any edge of the
// given face. With the restriction, that it can't be 'lastEdge'.
func calcHighestIntersection(v *Voronoi, bisector Edge, face FaceIndex, lastEdge EdgeIndex, lastVertex Vector) (EdgeIndex, Vector) {

    bestEdge := EmptyEdge
    bestIntersection := Vector{}

    minY := -100000000000.0
    maxY := -minY
    if lastVertex != (Vector{}) {
        maxY = lastVertex.Y
    }

    edge := v.faces[face].EEdge
    veryFirstEdge := edge

    first := true


    // Find be highest/best intersection of the bisector and edge.
    for edge != EmptyEdge && (first || edge != veryFirstEdge) {

        intersects, location := LineIntersection4(bisector, createLine(v, edge, true))
        //fmt.Printf("Found an intersection: %v, %v\n", intersects, location)

        if intersects &&
           edge != lastEdge &&
           location.Y < maxY &&
           location.Y > minY &&
           !Equal(location, lastVertex) &&
           !Equal(location, bestIntersection) {
            //fmt.Printf("--> Better than the last! %v\n", location)
            bestEdge = edge
            bestIntersection = location
            minY = bestIntersection.Y


        }

        first = false
        edge = v.edges[edge].ENext
    }

    return bestEdge, bestIntersection
}

// Calculates the upper or lower common support line for the two given convex hulls.
// Depending on the betterSide and betterPoint-function, it could potentially also calculate
// left/right common support lines, depending on the divide-and-conquer approach of the main algorithm.
func commonSupportLine(v *Voronoi, h1, h2 ConvexHull, betterSideL func(v1, v2 Vector) bool, betterSideR func(v1, v2 Vector) bool, betterPoint func(v1, v2, test Vector) bool) (FaceIndex, FaceIndex) {

    ai := h1.bestFaceIndex(v, betterSideL)
    bi := h2.bestFaceIndex(v, betterSideR)

    finished := false

    for !finished {
        finished = true

        nextAi := (ai-1)%len(h1)
        if nextAi < 0 {
            nextAi = len(h1)-1
        }

        // iterating through the convex hull points...
        if len(h1) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[nextAi]].ReferencePoint) {
            ai = nextAi
            finished = false
        }
        // iterating through the convex hull points...
        if len(h2) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[(bi+1)%len(h2)]].ReferencePoint) {
            bi = (bi+1)%len(h2)
            finished = false
        }

        // We try both directions every time! Because of some special cases, where the first point for h1 is the correct one, but for the very
        // first step (maybe later?), the next point is better for the moment. It then never goes back to the actual globally better point.
        // By checking both directions, we can actually get back. And there should be no downsides, as there is only ever one direction to go!

        // iterating through the convex hull points...
        if len(h1) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h1[(ai+1)%len(h1)]].ReferencePoint) {
            ai = (ai+1)%len(h1)
            finished = false
        }

        nextBi := (bi-1)%len(h2)
        if nextBi < 0 {
            nextBi = len(h2)-1
        }

        // iterating through the convex hull points...
        if len(h2) > 1 && betterPoint(v.faces[h1[ai]].ReferencePoint, v.faces[h2[bi]].ReferencePoint, v.faces[h2[nextBi]].ReferencePoint) {
            bi = nextBi
            finished = false
        }
    }

    return h1[ai], h2[bi]
}
// Helperfunctions for calculating the common support line.
// Were anonymous functions, but sourced out because they where getting too long.
func isBetterUp(v1, v2, test Vector) bool {
    // if the test-vector is left of the line, created by v1-v2.
    // Left means: on the upper side of the line (left).
    if IsLeft2D(v1, v2, test) {
        return true
    }
    side := SideOfLine(v1, v2, test)
    if math.Abs(side) <= EPS {
        return test.X >= v1.X && test.X <= v2.X
    }
    return false
}
func isBetterDown(v1, v2, test Vector) bool {
    if IsRight2D(v1, v2, test) {
        return true
    }

    side := SideOfLine(v1, v2, test)
    if math.Abs(side) <= EPS {
        //fmt.Printf("we get here.\n")
        return test.X > v1.X && test.X < v2.X
    }
    return false
}

//
// Extract the dividing chain between the two voronois without actually
// merging them! This will be done in an additional step.
//
func (v *Voronoi)extractDividingChain(left, right VoronoiEntryFace) []ChainElem {
    var dividingChain []ChainElem
    var loopCount = 0

    fmt.Printf("dividing chain for recursion: %v\n", g_recursions)

    h1 := v.ConvexHull(left)
    h2 := v.ConvexHull(right)

    fmt.Printf("h1: %v\n", h1)
    fmt.Printf("h2: %v\n", h2)

    isRight := func(v1, v2 Vector) bool {
        return v1.X > v2.X
    }
    isLeft := func(v1, v2 Vector) bool {
        return v1.X < v2.X
    }
    // Upper common support line!
    p, q := commonSupportLine(v, h1, h2, isLeft, isRight, isBetterUp)
    // Lower common support line!
    h1Down, h2Down := commonSupportLine(v,  h1, h2, isRight, isLeft, isBetterDown)

    // We don't cross the same edge twice!
    lastPEdge := EmptyEdge
    lastQEdge := EmptyEdge
    // last vertex we created on that separation line
    lastVertex := Vector{}

    bisector := PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
    //bisector = Amplify(bisector, 100.0)

    // As long as we didn't reach the lowest possible tangente, we continue.
    for {
        if g_recursions == 156 && loopCount == 160 {
            //fmt.Printf("draw dividing chain. LoopCount == 200\n")
            drawDividingChain(dividingChain)
            os.Exit(0)
        }
        //fmt.Printf("for with p == %v, q == %v\n", p, q)
        // Here we break out of the loop, when we reach the very bottom!
        lastMerge := p == h1Down && q == h2Down
        fmt.Printf("p: %v, q: %v, h1Down: %v, h2Down: %v\n", p, q, h1Down, h2Down)

        edgeP, locationP := calcHighestIntersection(v, bisector, p, lastPEdge, lastVertex)
        fmt.Printf("found P intersection: %v\n", locationP.Y)
        edgeQ, locationQ := calcHighestIntersection(v, bisector, q, lastQEdge, lastVertex)
        fmt.Printf("found Q intersection: %v\n", locationQ.Y)

        switch {

            // For the case, that we merge two trivial voronois with no edges or
            // the very last step. Now just create two edges and we're done.
            case edgeP == EmptyEdge && edgeQ == EmptyEdge:
                fmt.Printf("e\n")
                dividingChain = append(dividingChain, ChainElem{Vector{}, EmptyEdge, EmptyEdge, p, q, bisector})

            // Equal Intersection with both edges. So 4 edges meet.
            case edgeQ != EmptyEdge && edgeP != EmptyEdge && Equal(locationP, locationQ):
                fmt.Printf("p && q\n")
                // This could result in kind of an ambiguous Delaunay triangulations. At least regarding a triangle.
                // This should now work for situations, where 4 edges meet. I have to re-examine for even more edges meeting...

                dividingChain = append(dividingChain, ChainElem{locationP, edgeP, edgeQ, p, q, bisector})

                lastVertex   = locationP
                lastPEdge    = edgeP
                lastQEdge    = edgeQ

                p = v.edges[v.edges[edgeP].ETwin].FFace
                q = v.edges[v.edges[edgeQ].ETwin].FFace

            // We intersect with an edge of the face p
            case edgeP != EmptyEdge && (edgeQ == EmptyEdge || locationP.Y >= locationQ.Y):
                //fmt.Printf("p\n")
                dividingChain = append(dividingChain, ChainElem{locationP, edgeP, EmptyEdge, p, q, bisector})

                lastVertex   = locationP
                lastPEdge    = edgeP
                lastQEdge    = EmptyEdge

                p = v.edges[v.edges[edgeP].ETwin].FFace

            // We intersect with an edge of the face q
            case edgeQ != EmptyEdge && (edgeP == EmptyEdge || locationQ.Y >= locationP.Y):
                //fmt.Printf("q\n")
                dividingChain = append(dividingChain, ChainElem{locationQ, EmptyEdge, edgeQ, p, q, bisector})

                lastVertex   = locationQ
                lastQEdge    = edgeQ
                lastPEdge    = EmptyEdge

                q = v.edges[v.edges[edgeQ].ETwin].FFace


            default:
                fmt.Println("default")
                lastMerge = true
                v.pprint()
                os.Exit(0)

        }


        if lastMerge {
            break
        }

        bisector = PerpendicularBisector(v.faces[p].ReferencePoint, v.faces[q].ReferencePoint)
        //bisector = Amplify(bisector, 100.0)
        loopCount += 1
    }

    if g_recursions == 24 {
        //dividingChain(dividingChain)
    }

    return dividingChain
}

// Right now, a voronoi diagram is identified by a Vertex.
// Here two not overlapping voronoi diagrams are merged.
// They HAVE to be left/right of each other with NO overlapping. This HAS to be guaranteed!
func (v *Voronoi)mergeVoronoi(left, right VoronoiEntryFace) VoronoiEntryFace {

    g_recursions += 1;

    nextPEdge    := EmptyEdge
    lastDownEdge := EmptyEdge
    lastVertex   := EmptyVertex

    lastQEdge    := EmptyEdge
    lastPEdge    := EmptyEdge

    chainList := v.extractDividingChain(left, right)

    //fmt.Printf("chainList: %v\n", chainList)

    // Iterating through the dividing chain and actually merge the voronois.
    for i,chain := range chainList {

        //fmt.Printf("i == %v\n", i)
        heVertex := EmptyVertex
        if chain.intersection != (Vector{}) {

            if chain.intersection == InfinitePoint {
                fmt.Println("INFINITY")
                heVertex = InfiniteVertex
            } else {
                heVertex = v.createVertex(chain.intersection)
            }
        }

        otherWayBisector := chain.bisector
        otherWayBisector.Dir = Mult(otherWayBisector.Dir, -1)

        heEdgeUp   := v.createEdge(heVertex,   EmptyEdge, chain.edgeP,  nextPEdge,   chain.p, otherWayBisector)
        heEdgeDown := v.createEdge(lastVertex, heEdgeUp,  lastDownEdge, chain.edgeQ, chain.q, chain.bisector)
        v.edges[heEdgeUp].ETwin = heEdgeDown

        //fmt.Printf("heEdgeUp  : %v\n", v.edges[heEdgeUp])
        //fmt.Printf("heEdgeDown: %v\n", v.edges[heEdgeDown])

        if lastDownEdge != EmptyEdge {
            v.edges[lastDownEdge].ENext = heEdgeDown
        }

        // Set Previous edge. Only in case of an intersection with Q. The P one is handled implicitely.
        if nextPEdge != EmptyEdge {
            v.edges[nextPEdge].EPrev = heEdgeUp
        }

        // For the very first dividing chain element holds, that the face q has to update
        // its edge reference to this down-edge. No matter what.
        if i == 0 {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if i == len(chainList)-1 {
            v.faces[chain.p].EEdge = heEdgeUp
        }

        // If edgeQ is the current edge for the face q and we have an intersection with
        // edgeQ, the new first edge for the face q must be edgeQ!
        // Simple case!
        if chain.edgeQ != EmptyEdge && v.faces[chain.q].EEdge == chain.edgeQ && chain.intersection != InfinitePoint {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if v.faces[chain.q].EEdge == EmptyEdge || (lastDownEdge == EmptyEdge && heVertex != InfiniteVertex) {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if v.faces[chain.q].EEdge == chain.edgeQ && chain.intersection != InfinitePoint {
            v.faces[chain.q].EEdge = heEdgeDown
        }

        if v.faces[chain.p].EEdge == EmptyEdge || !heVertex.Valid() {
            v.faces[chain.p].EEdge = heEdgeUp
        }

        // Removing stuff!
        if chain.edgeP != EmptyEdge {
            // Delete vertices that overlap with the other side (to the right)!
            vertex := v.edges[v.edges[chain.edgeP].ETwin].VOrigin

            if vertex.Valid() && v.vertices[vertex] != emptyV {
                //fmt.Printf("    ----> Found a vertex that has do be deleted (P) (v-index: %v - %v)\n", vertex, v.vertices[vertex].Pos)

                removeNextEdge := v.edges[chain.edgeP].ENext

                // Now, this is somewhat tricky. We say, we want to keep an O(n logn) runtime in general.
                // The divide&conquer takes O(logn), while the merging of two voronois can take up O(n) time
                // as we have to go through every face, if the dividing chain goes through all of them.
                // Now this loop here can run a maximum of O(n) as well, which logically can add up to O(n^2).
                // Fortunately, this is not the case, because this loop here can go through O(n) edges only over
                // one complete merge. It also deletes all n objects, it touches; meaning, that if we remove edges,
                // they are not relevant for the merging any more. So both loops can only add up to one O(n).
                for removeNextEdge != EmptyEdge && v.edges[removeNextEdge] != emptyE && v.edges[removeNextEdge].ETwin != lastPEdge {
                    tmpRemoveNextEdge := v.edges[removeNextEdge].ENext

                    tmpF := v.edges[v.edges[removeNextEdge].ETwin].FFace
                    if v.faces[tmpF].EEdge == v.edges[removeNextEdge].ETwin {
                        v.faces[tmpF].EEdge = v.edges[v.edges[removeNextEdge].ETwin].ENext
                    }

                    if v.faces[v.edges[removeNextEdge].FFace].EEdge == removeNextEdge {
                        v.faces[v.edges[removeNextEdge].FFace].EEdge = chain.edgeP
                    }

                    fmt.Printf("DELETE E: %v | %v\n", removeNextEdge, v.edges[removeNextEdge].ETwin)

                    v.deleteVertex(v.edges[removeNextEdge].VOrigin)
                    v.deleteEdgePair(v.edges[removeNextEdge].ETwin, removeNextEdge)

                    removeNextEdge = tmpRemoveNextEdge

                }
            }
        }

        // Removing Stuff
        if chain.edgeQ != EmptyEdge {
            // Delete vertices that overlap with the other side (to the left)!
            vertex := v.edges[chain.edgeQ].VOrigin
            if vertex.Valid() && v.vertices[vertex] != emptyV {
                //fmt.Printf("    ----> Found a vertex that has do be deleted (Q) (v-index: %v)\n", vertex)

                removeLastEdge := v.edges[chain.edgeQ].EPrev

                for removeLastEdge != EmptyEdge && v.edges[removeLastEdge] != emptyE && v.edges[removeLastEdge].ETwin != lastQEdge {

                    tmpRemoveLastEdge := v.edges[removeLastEdge].EPrev

                    tmpF := v.edges[v.edges[removeLastEdge].ETwin].FFace
                    if v.faces[tmpF].EEdge == v.edges[removeLastEdge].ETwin {
                        v.faces[tmpF].EEdge = v.edges[v.edges[removeLastEdge].ETwin].EPrev
                    }

                    if v.faces[v.edges[removeLastEdge].FFace].EEdge == removeLastEdge {
                        v.faces[v.edges[removeLastEdge].FFace].EEdge = heEdgeDown
                    }

                    v.deleteVertex(v.edges[v.edges[removeLastEdge].ETwin].VOrigin)
                    v.deleteEdgePair(v.edges[removeLastEdge].ETwin, removeLastEdge)

                    removeLastEdge = tmpRemoveLastEdge

                }

            }
        }

        // Setting States
        if chain.edgeP != EmptyEdge {
            lastPEdge = chain.edgeP

            v.edges[chain.edgeP].ENext = heEdgeUp
            v.edges[v.edges[chain.edgeP].ETwin].VOrigin = heVertex
            nextPEdge = v.edges[chain.edgeP].ETwin

        } else {
            nextPEdge = heEdgeUp
        }

        // Setting States
        if chain.edgeQ != EmptyEdge {
            v.edges[chain.edgeQ].EPrev = heEdgeDown
            lastQEdge = chain.edgeQ
            v.edges[chain.edgeQ].VOrigin = heVertex
            lastDownEdge = v.edges[chain.edgeQ].ETwin

        } else {
            lastDownEdge = heEdgeDown
        }

        lastVertex = heVertex

    }

    fmt.Printf("FINISHED MERGE OF %v AND %v\n\n", left, right)

    return left
}

// Voronoi divide and conquer entry point
func (v *Voronoi)divideAndConquer(points PointList) VoronoiEntryFace {
    l := len(points)

    // Recursion stops at one point!
    if l == 1 {
        // One face without anything else is the most primitive, trivial Voronoi tessellation!
        return VoronoiEntryFace(v.createFace(points[0], EmptyEdge))
    }

    return v.mergeVoronoi(v.divideAndConquer(points[:l/2]), v.divideAndConquer(points[l/2:]))
}

// Sorts the points and returns a voronoi tessellation.
func CreateVoronoi(pointList PointList) Voronoi {
    sort.Sort(pointList)
    n := len(pointList)

    // See: http://www.cs.wustl.edu/~pless/546/lectures/L11.html
    // for the calculations of maximum voronoi object count.
    v := Voronoi {
        vertices:           make([]HEVertex, 2*n-5 + 3),
        firstFreeVertexPos: 0,
        edges:              make([]HEEdge, 2*(3*n-6) + 6),
        firstFreeEdgePos:   0,
        faces:              make([]HEFace, n),
        firstFreeFacePos:   0,
    }

    //for i,_ := range v.vertices {
    //    v.vertices[i] = emptyV
    //}
    //for i,_ := range v.edges {
    //    v.edges[i] = emptyE
    //}
    //for i,_ := range v.faces {
    //    v.faces[i] = emptyF
    //}

    g_freeEdgePositions.Init()
    g_freeVertexPositions.Init()

    v.divideAndConquer(pointList)

    e := v.Verify()

    if e != nil {
        fmt.Printf("VERIFICATION_ERROR: '%v'\n", e)
    }

    return v
}

//
// Normal test distribution of points. No special cases for Voronoi generation.
//
func testNormal01() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 50., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_01", true)

}

func testNormal02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 30., 0})
    pointList = append(pointList, Vector{80., 20., 0})
    pointList = append(pointList, Vector{55., 40., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_02", true)
}

func testNormal03() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_normal_03\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{45., 45., 0})
    pointList = append(pointList, Vector{50., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_normal_03", true)

}

//
// Test cases where 4 lines meet and p/q have the same intersection point.
//
func testEqualIntersection01() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_intersection_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{28, 30, 0})
    pointList = append(pointList, Vector{30, 20, 0})
    pointList = append(pointList, Vector{30, 40, 0})

    pointList = append(pointList, Vector{40, 20, 0})
    pointList = append(pointList, Vector{40, 40, 0})
    pointList = append(pointList, Vector{42, 30, 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_equal_intersection_01", true)
}

func testEqualIntersection02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_intersection_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 50., 0})
    pointList = append(pointList, Vector{40., 20., 0})
    pointList = append(pointList, Vector{60., 20., 0})
    pointList = append(pointList, Vector{70., 50., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_equal_intersection_02", true)
}

//
// Test cases with linear dependent points or lines
//
func testLinearDepentence0() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_0\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20, 70, 0})
    pointList = append(pointList, Vector{30, 70, 0})
    pointList = append(pointList, Vector{40, 70, 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_01a", true)
}

//
// Test cases with linear dependent points or lines
//
func testLinearDepentence01a() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20, 70, 0})
    pointList = append(pointList, Vector{30, 70, 0})
    pointList = append(pointList, Vector{40, 70, 0})
    pointList = append(pointList, Vector{50, 70, 0})

    //pointList[0].Y -= 1;
    //pointList[3].Y -= 1;

    //pointList = pointList[:len(pointList)/2]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_01a", true)
}

//
// Special case, merging two linear dependent voronoi with eath other!
//
func testLinearDepentence05() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_05\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20, 70, 0})
    pointList = append(pointList, Vector{30, 70, 0})
    pointList = append(pointList, Vector{40, 70, 0})

    pointList = append(pointList, Vector{20, 80, 0})
    pointList = append(pointList, Vector{30, 80, 0})
    pointList = append(pointList, Vector{40, 80, 0})
    //pointList = append(pointList, Vector{50, 70, 0})

    //pointList[0].Y -= 1;
    //pointList[3].Y -= 1;

    //pointList = pointList[:len(pointList)/2]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_01a", true)
}

func testLinearDepentence01b() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_01\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20, 70, 0})
    pointList = append(pointList, Vector{30, 70, 0})
    pointList = append(pointList, Vector{40, 70, 0})
    pointList = append(pointList, Vector{50, 70, 0})

    //pointList[0].Y -= 1;
    //pointList[3].Y -= 1;

    //pointList = pointList[:len(pointList)/2]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_01b", true)
}

func testLinearDepentence02() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_02\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_02", true)
}

func testLinearDepentence03() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_03\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{20., 10., 0})
    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})
    pointList = append(pointList, Vector{50., 10., 0})
    pointList = append(pointList, Vector{60., 10., 0})
    pointList = append(pointList, Vector{70., 10., 0})
    pointList = append(pointList, Vector{80., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_03", true)
}

func testLinearDepentence04() {
    fmt.Printf("=============================\n")
    fmt.Printf("=== test_linear_dependence_04\n")
    fmt.Printf("=============================\n")
    var pointList PointList

    pointList = append(pointList, Vector{30., 10., 0})
    pointList = append(pointList, Vector{40., 10., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_linear_dependence_04", true)
}

//
// Test cases that cause infinite loops or unknown problems.
//

func testUnknownProblem01() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_01\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 15., 0})
    pointList = append(pointList, Vector{15., 65., 0})
    pointList = append(pointList, Vector{20., 40., 0})

    pointList = append(pointList, Vector{40.,  5., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 60., 0})
    // Adding this point results in a wrong voronoi! Really messed up.
    pointList = append(pointList, Vector{90., 10., 0})

    for i,_ := range pointList {
        pointList[i] = Mult(pointList[i], 0.5)
        pointList[i] = Add(pointList[i], Vector{25,25,0})
    }

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_01", true)
}

func testUnknownProblem02() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_02\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 10., 0})
    pointList = append(pointList, Vector{20., 20., 0})

    pointList = append(pointList, Vector{25., 40., 0})

    pointList = append(pointList, Vector{30., 10., 0})

    pointList = append(pointList, Vector{40., 30., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 10., 0})


    //pointList = pointList[:len(pointList)/2]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_02", true)
}

func testUnknownProblem03() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_03\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{10., 15., 0})
    pointList = append(pointList, Vector{15., 65., 0})
    pointList = append(pointList, Vector{20., 40., 0})

    pointList = append(pointList, Vector{40.,  5., 0})
    pointList = append(pointList, Vector{50., 20., 0})
    pointList = append(pointList, Vector{60., 60., 0})

    // Infinite loop by adding this point.
    pointList = append(pointList, Vector{70., 70., 0})

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_03", true)
}

func testUnknownProblem04() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_04\n")
    fmt.Printf("===========================\n")
    var pointList PointList

    pointList = append(pointList, Vector{26.140400405282577, 61.621512929924580, 0})
    pointList = append(pointList, Vector{30.777003166789225, 25.899889187518045, 0})
    pointList = append(pointList, Vector{38.049244986921730, 50.778804685015940, 0})
    pointList = append(pointList, Vector{39.182654173412490, 53.260689369108890, 0})
    pointList = append(pointList, Vector{56.653551080486500, 67.211346445470530, 0})

    pointList = append(pointList, Vector{61.682906469053240, 71.811352995110340, 0})
    pointList = append(pointList, Vector{61.787090244227150, 57.606346468267920, 0})
    pointList = append(pointList, Vector{69.598081344490680, 56.223764727531960, 0})
    pointList = append(pointList, Vector{71.996875920571910, 27.832275576849792, 0})
    pointList = append(pointList, Vector{72.329734102573130, 63.247795405728404, 0})


    if false {
        for i,_ := range pointList {
            //pointList[i] = Mult(pointList[i], 0.5)
            pointList[i] = Add(pointList[i], Vector{-20,0,0})
        }
    }

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_04", true)
}

func testUnknownProblem05() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_05\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409032579303917
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_05", true)
}

func testUnknownProblem06() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_06\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409509781449406
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }
    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_06", true)
}

func testUnknownProblem07() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_07\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409694089411269
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }
    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_07", true)
}

func testUnknownProblem08() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_08\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409796721354066
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_08", true)
}

func testUnknownProblem09() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_09\n")
    fmt.Printf("===========================\n")
    count := 5
    var seed int64 = 1482476817056249360
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        //v := Vector{r.Float64()*150., r.Float64()*150., 0}
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[2:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_09", false)
}

func testUnknownProblem10() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_10\n")
    fmt.Printf("===========================\n")
    count := 5
    var seed int64 = 1482477466512747360
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)

    //pointList = pointList[:5]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_10", true)
}

func testUnknownProblem11() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_11\n")
    fmt.Printf("===========================\n")
    count := 10
    var seed int64 = 1482409032579303917
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)
    pointList = pointList[:8]

    //pointList = pointList[:4]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_11", true)
}

func testUnknownProblem12() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_12\n")
    fmt.Printf("===========================\n")
    count := 9
    var seed int64 = 1483370038194119290
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*50.+25., r.Float64()*50.+25., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)
    pointList = pointList[1:]

    //pointList = pointList[4:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_12", true)
}

func testUnknownProblem13() {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_unknown_problem_13\n")
    fmt.Printf("===========================\n")
    count := 10000
    var seed int64 = 1499748501205604561
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*100., r.Float64()*100., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)
    //pointList = pointList[1:]

    pointList = pointList[:len(pointList)/2]
    pointList = pointList[:len(pointList)/2]
    pointList = pointList[len(pointList)/2:]
    pointList = pointList[len(pointList)/2:]
    pointList = pointList[len(pointList)/2:]
    pointList = pointList[len(pointList)/2:]

    pointList = pointList[len(pointList)/2:]
    pointList = pointList[len(pointList)/2:]

    // Before coming here, both *100 and *500 are perfectly identical!
    // Not true, it just doesn't cause an error...

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_13", true)
}

func testUnknownProblemSeed(seed int64, count int) {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_seed_%v\n", strconv.FormatInt(seed, 10))
    fmt.Printf("===========================\n")
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*1000., r.Float64()*1000., 0}
        pointList = append(pointList, v)
    }

    sort.Sort(pointList)


    for _,v := range pointList {
        fmt.Printf("%v\n", v)
    }

    //pointList = pointList[len(pointList)/2:]
    //pointList = pointList[len(pointList)/2:]

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_unknown_problem_seed_" + strconv.FormatInt(seed, 10), true)
}

//
// Test cases with random points
//

func testRandom(count int) {
    fmt.Printf("===========================\n")
    fmt.Printf("=== test_random_%v\n", count)
    fmt.Printf("===========================\n")
    seed := time.Now().UTC().UnixNano()
    r := rand.New(rand.NewSource(seed))
    var pointList PointList

    fmt.Printf("Seed: %v\n", strconv.FormatInt(seed, 10))
    fmt.Printf("Random poins:\n\t")
    for i:= 0; i < count; i++ {
        v := Vector{r.Float64()*100., r.Float64()*100., 0}
        pointList = append(pointList, v)
        fmt.Printf("%v, ", v)
    }
    fmt.Printf("\n\n")

    v := CreateVoronoi(pointList)
    v.pprint()

    ch := v.ConvexHull(0)
    fmt.Println(ch)

    v.createImage("test_random_" + strconv.FormatInt(seed, 10), true)
}


func main() {

    g_drawImages = false
    g_pprint = false
    working := false

    if working {
        testNormal01()
        testNormal02()
        testNormal03()
        testEqualIntersection01()
        testEqualIntersection02()
        testLinearDepentence04()
        testUnknownProblem03()
        testUnknownProblem09()
        testUnknownProblem10()
        testUnknownProblemSeed(1489941049442429888, 5)
        testUnknownProblem01()
        testUnknownProblem08()
        testUnknownProblem11()
        testUnknownProblem06()
        testUnknownProblem07()
        testUnknownProblem12()
        testUnknownProblem04()
        testUnknownProblemSeed(1483369884537650258, 20)
        testUnknownProblemSeed(1483370089898481236, 15)
        testUnknownProblem05()
        testUnknownProblem02()
        testUnknownProblemSeed(1483370150842201370, 15)
        testUnknownProblemSeed(1483370130545841965, 15)
        testUnknownProblemSeed(1496121738043120503, 20)
        testUnknownProblemSeed(1496247478836154924, 40)
        testUnknownProblemSeed(1496247359069208740, 50)
        testUnknownProblemSeed(1496294040279258517, 100)
        testUnknownProblemSeed(1496294255661064215, 500)
        testUnknownProblemSeed(1496295030863743777, 10000)
    }

    toBeVerified := false

    if toBeVerified {

        testLinearDepentence01a()
        // works.
        // A first edge {-1 0 6 1 {{35 510 0} {0 -1000 0}}}: 1 must be referenced as first edge by the corresponding face 1 !
        fmt.Println("Test: testLinearDepentence02")
        testLinearDepentence0()
        testLinearDepentence02()
        testLinearDepentence05()

        // works.
        // A first edge {-2 1 -1 1 {{35 510 0} {-0 1000 -0}}}: 0 must be referenced as first edge by the corresponding face 1 !
        fmt.Println("Test: testLinearDepentence03")
        testLinearDepentence03()

    }

    test := true
    g_pprint = true
    g_drawImages = true

    if test {
        //testRandom(10000)

        //testUnknownProblemSeed(1499748501205604561, 10000)
        testUnknownProblem13()

        //testLinearDepentence0()
        //testLinearDepentence01b()
    }

}




