package vector

import (
    "fmt"
)

// Because of laziness, a Vector could also just be a point. Depending on context.
type Vector struct {
    X,Y,Z   float32
}

type Edge struct {
    Pos             Vector
    Dir             Vector
}

// Implement the sort.Interface, so it can be sorted by the standard algorithm of Go.
type PointList []Vector
func (slice PointList) Len() int {
    return len(slice)
}
func (slice PointList) Less(i, j int) bool {
    return slice[i].X < slice[j].X;
}
func (slice PointList) Swap(i, j int) {
    slice[i], slice[j] = slice[j], slice[i]
}

// Just makes the edge a lot bigger according to s.
func (e Edge) Amplify(s float32) {
    e.Dir.Mult(2.0*s)
    e.Pos = e.Pos.Sub(e.Dir)
}

func (e Edge) Copy() Edge {
    return Edge{Pos: e.Pos.Copy(), Dir: e.Dir.Copy()}
}


func Add(v1, v2 Vector) Vector {
    return Vector{X: v1.X+v2.X, Y: v1.Y+v2.Y, Z: v1.Z+v2.Z}
}
func Mult(v1 Vector, s float32) Vector {
    return Vector{X: v1.X*s, Y: v1.Y*s, Z: v1.Z*s}
}
func Sub(v1, v2 Vector) Vector {
    return Vector{X: v2.X-v1.X, Y: v2.Y-v1.Y, Z: v2.Z-v1.Z}
}

func (v Vector) Copy() Vector {
    return Vector{X:v.X, Y:v.Y, Z:v.Z}
}

func (v Vector) Add(v1 Vector) Vector {
    v.X += v1.X
    v.Y += v1.Y
    v.Z += v1.Z
    return v
}

func (v Vector) Sub(v1 Vector) Vector {
    v.X -= v1.X
    v.Y -= v1.Y
    v.Z -= v1.Z
    return v
}

func (v Vector) Mult(s float32) Vector {
    v.X *= s
    v.Y *= s
    v.Z *= s
    return v
}

// Calculates v2 - v1.
func Diff (v1, v2 Vector) Vector {
    return Vector {
        X: v2.X - v1.X,
        Y: v2.Y - v1.Y,
        Z: v2.Z - v1.Z,
    }
}

// Calculates the cross product: v1 x v2
func Cross(v1, v2 Vector) Vector {
    return Vector {
        X: v1.Y*v2.Z - v1.Z*v2.Y,
        Y: v1.Z*v2.X - v1.X*v2.Z,
        Z: v1.X*v2.Y - v1.Y*v2.X,
    }
}

// Calculates the middle point between p1 and p2.
func MiddlePoint(p1, p2 Vector) Vector {
    d := Diff(p1, p2)
    return Vector {
        X: p1.X + d.X/2.0,
        Y: p1.Y + d.Y/2.0,
        Z: p1.Z + d.Z/2.0,
    }
}

// Calculates an edge representation of one perpendicular bisector between p1 and p2
func PerpendicularBisector(p1, p2 Vector) Edge {
    return Edge {
        Pos: MiddlePoint(p1, p2),
        Dir: Cross(Diff(p1, p2), Vector{0,0,1}),
    }
}

// Determinante of v1 and v2, ignoring the Z-component!
// Honestly, I don't really know, what this does^^
func Det2D(v1, v2 Vector) float32 {
    return v1.X*v2.Y - v1.Y*v2.X
}

// Calculates the intersection point of v1 and v2, if it exists.
// http://stackoverflow.com/questions/20677795/find-the-point-of-intersecting-lines
func LineIntersection(e1 Edge, e2 Edge) Vector {

    xdiff := Vector{X: e1.Dir.X, Y: e2.Dir.X}
    ydiff := Vector{X: e1.Dir.Y, Y: e2.Dir.Y}

    div := Det2D(xdiff, ydiff)
    if div <= 0.0001 {
        fmt.Println("Lines do NOT intersect!")
        return Vector{}
    }

    d := Vector{X: Det2D(e1.Pos, Add(e1.Pos, e1.Dir)), Y: Det2D(e2.Pos, Add(e2.Pos, e2.Dir))}

    return Vector{
        X: Det2D(d, xdiff) / div,
        Y: Det2D(d, ydiff) / div,
        Z: 0.0,
    }
}



