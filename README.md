# Approximate 2d point set by k vertex convex polygon

### Written By Sariel Har-Peled, 2025-Sep-16


This code computes a polygon that approximates the convex-hull of a
given set of points (all the in the plane). Specifically, one is given
a parameter k, and the algorithm adds the k points that each of them
add the most area to the current convex-hull. For simplicity, the two
starting points are the two x-extreme points. The algorithms works in
O(n log n) time if you are lucky, and O(n (log n + k) ) if you are
not. There are ways to improve the running time to O(n log n) but I
doubt if it is worth the effort.

Note that this is a heuristic. Computing the k-polygon maximizing the
area can be done by dynamic programming but is significantly more
tedious and slower (i.e., O(n^3)?). This hack should be good enough in
practice.

Most of the code was written by prompting Gemini to generate some
pieces of code and putting them together (with a bit of coding myself
[ha, the suffering]).

