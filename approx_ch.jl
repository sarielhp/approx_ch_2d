#! /bin/env julia

# Written By Sariel Har-Peled, 2025-Sep-16
#
# This code computes a polygon that approximates the convex-hull of a
# given set of points (all the in the plane). Specifically, one is
# given a parameter k, and the algorithm adds the k points that each
# of them add the most area to the current convex-hull. For
# simplicity, the two starting points are the two x-extreme
# points. The algorithms works in O(n log n) time if you are lucky,
# and O(n (log n + k) ) if you are not. There are ways to improve the
# running time to O(n log n) but I doubt if it is worth the effort.
#
# Note that this is a heuristic. Computing the k-polygon maximizing
# the area can be done by dynamic programming but is significantly
# more tedious and slower (i.e., O(n^3)?). This hack should be good
# enough in practice.
#
# Most of the code was written by prompting Gemini to generate
# some pieces of code and putting them together (with a bit of coding
# myself [ha, the suffering]).

using Random, Cairo, DataStructures

struct Point
    x::Float64
    y::Float64
end


"""
    is_left_turn(p1::Point, p2::Point, p3::Point)

Checks if the sequence of points p1, p2, p3 forms a left turn.
A left turn is indicated by a positive signed area of the triangle formed by the three points.
This is also known as a counter-clockwise turn.
"""
function is_left_turn(p1::Point, p2::Point, p3::Point)::Bool
    # The signed area formula is 0.5 * (x1(y2-y3) + x2(y3-y1) + x3(y1-y2)).
    # We can omit the 0.5 since we only care about the sign.
    cross_product = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)

    # A positive cross product indicates a counter-clockwise (left) turn.
    return cross_product > 0
end
function is_leq_turn(p1::Point, p2::Point, p3::Point)::Bool
    # The signed area formula is 0.5 * (x1(y2-y3) + x2(y3-y1) + x3(y1-y2)).
    # We can omit the 0.5 since we only care about the sign.
    cross_product = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)

    # A positive cross product indicates a counter-clockwise (left) turn.
    return cross_product >= 0.0
end


function convex_hull( _points::Vector{Point} )
    n = length( _points )
    if n <= 3
        return unique(_points) # Handle small cases and duplicates
    end

    # Step 1: Sort points by x-coordinate, with y as a tie-breaker
    points = sort(_points, by = p -> (p.x, p.y))

    # Step 2: Build the lower hull
    lower_hull = Point[]
    for p in points
        while length(lower_hull) >= 2 && is_leq_turn(lower_hull[end-1], lower_hull[end], p)
            pop!(lower_hull)
        end
        push!(lower_hull, p)
    end

    # Step 3: Build the upper hull
    upper_hull = Point[]
    for i in n:-1:1
        p = points[i]
        while length(upper_hull) >= 2 && is_leq_turn(upper_hull[end-1], upper_hull[end], p)
            pop!(upper_hull)
        end
        push!(upper_hull, p)
    end

    # Step 4: Combine the two hulls, removing duplicates at the start and end
    # The last point of the lower hull and the last point of the upper hull are the same
    # as the first point of the upper and lower hulls, respectively.
    pop!(lower_hull)
    pop!(upper_hull)
    
    return [lower_hull; upper_hull]
end


function sort_points_by_x(points::Vector{Point})
    return sort(points, by = p -> p.x)
end



function generate_random_points_in_disk(n::Int, radius::Float64)
    points = Vector{Point}(undef, n)
    for i in 1:n
        angle = 2π * rand()
        r = radius * sqrt(rand())

        # Convert polar coordinates (r, angle) to Cartesian coordinates (x, y)
        x = r * cos(angle)
        y = r * sin(angle)
        points[i] = Point( x, y )
    end
    return points
end

function generate_random_points(n::Int; x_range=(0.0, 1.0), y_range=(0.0, 1.0))
    points = Vector{Point}(undef, n)

    for i ∈ 1:n
        rand_x = rand() * (x_range[2] - x_range[1]) + x_range[1]
        rand_y = rand() * (y_range[2] - y_range[1]) + y_range[1]
        
        points[i] = Point(rand_x, rand_y)
    end

    return points
end


struct  SubProb
    i::Int64
    middle::Int64
    j::Int64    
    val::Float64
end

function triangle_area(p1::Point, p2::Point, p3::Point)
    # The shoelace formula for a triangle:
    # Area = 0.5 * |(x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2))|
    # This is equivalent to half the magnitude of the cross product of two vectors
    # formed by the sides of the triangle.
    area = 0.5 * abs((p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)))
    return area
end


function   push_candidate_sol( heap, P, i, m, j )
    if  ( m == 0 )
        return
    end
        
    area = triangle_area( P[ i ], P[m], P[ j ] )
    enqueue!( heap, SubProb( i, m, j, area), area )
end


function    get_best_point( P, i, j )
    if  ( abs( i - j ) <= 1 )
        return  0
    end

    step_val = (i <= j) ? 1 : -1
    val = -1.0
    i_best = 0
    for  x ∈ i+step_val:step_val:j-step_val
        if  ! is_left_turn( P[i], P[x], P[j] )
            continue
        end
        area = triangle_area( P[i], P[x], P[j] )
        if  area > val
            val = area
            i_best = x
        end
    end
    return  i_best
end


function  compute_good_approx( _P::Vector{Point}, k::Int )
    n = length( _P )
    if  ( n <= max( k, 3 ) )
        return  convex_hull( _P )
    end

    P = sort_points_by_x( _P )
    
    out = Vector{Point}()
    heap = PriorityQueue{SubProb, Float64}(Base.Order.Reverse)
    
    push!( out, first( P ) );
    push!( out, last( P ) );
    k -= 2

    i_a = get_best_point( P, 1, length( P ) )
    i_b = get_best_point( P, length( P ), 1 )

    push_candidate_sol( heap, P, 1, i_a, length( P ) )
    push_candidate_sol( heap, P, length( P ), i_b, 1 )
    while  ( ! isempty( heap ) )  &&  ( k > 0 ) 
        sub, val = peek( heap )
        dequeue!( heap )

        push!( out, P[ sub.middle ] );
        k -= 1

        i_m = get_best_point( P, sub.i, sub.middle )
        if  ( i_m > 0 )
            push_candidate_sol( heap, P, sub.i, i_m, sub.middle )
        end

        j_m = get_best_point( P, sub.middle, sub.j )
        if  ( i_m > 0 )
            push_candidate_sol( heap, P, sub.middle, j_m, sub.j )
        end
    end

    return  convex_hull( out )
end


function draw_polygon(vertices::Vector{Point}, ctx )
    Cairo.move_to(ctx, vertices[1].x, vertices[1].y)

    # Draw a line to each subsequent vertex
    for i in 2:length(vertices)
        Cairo.line_to(ctx, vertices[i].x, vertices[i].y)
    end

    # Close the path by connecting the last vertex back to the first
    Cairo.close_path(ctx)

    # Set the color for the fill and stroke
    Cairo.set_source_rgb(ctx, 0.2, 0.4, 0.8) # A nice blue color for filling
    
    # Fill the polygon with the current color
    Cairo.fill_preserve(ctx) # `fill_preserve` keeps the path for stroking

    # Set the color and width for the outline
    Cairo.set_source_rgb(ctx, 0.0, 0.0, 0.0) # Black for the outline
    Cairo.set_line_width(ctx, 0.01)

    # Stroke (draw) the outline of the polygon
    Cairo.stroke(ctx)
end

function draw_points_polygon_to_pdf(points::Vector{Point}, out, filename::String,
                                    disk_radius::Float64)
    # Define PDF page dimensions
    width = 600
    height = 600

    # Create a PDF drawing surface
    surface = Cairo.CairoPDFSurface(filename, width, height)
    ctx = Cairo.CairoContext(surface)

    # Set up a transformation to center the drawing and scale it
    # We add a small margin (e.g., 10%)
    margin = 1.1
    Cairo.set_line_width(ctx, 1.0)
    Cairo.set_source_rgb(ctx, 0, 0, 0) # Black
    Cairo.translate(ctx, width / 2, height / 2)
    Cairo.scale(ctx, width / (2 * disk_radius * margin), -height / (2 * disk_radius * margin))

    draw_polygon( out, ctx )

    
    # Draw the boundary of the disk
    #Cairo.set_line_width(ctx, 0.01)
    #Cairo.arc(ctx, 0, 0, disk_radius, 0, 2π)
    #Cairo.stroke(ctx)

    # Draw each point
    Cairo.set_source_rgb(ctx, 0.7, 0.2, 0.2) # A red color for the points
    for p in points
        # Draw a small circle for each point
        Cairo.arc(ctx, p.x, p.y, disk_radius * 0.01, 0, 2π) # Size of the point
        Cairo.fill(ctx)
    end

    # Clean up and save the PDF
    Cairo.finish(surface)
end


function  (@main)(args)
    n = 400;
    k = 7;
    radius = 10.0;
    
    P = generate_random_points_in_disk( n, radius )
    out = compute_good_approx( P, k )
    @assert( length( out ) == k )
    
    draw_points_polygon_to_pdf( P, out, "test.pdf", radius )
    
    return  0
end
