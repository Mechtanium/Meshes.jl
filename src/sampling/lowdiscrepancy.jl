
# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#TODO: Update the docstring for LowDiscrepancySampling below
"""
    HomogeneousSampling(size, [weights])

Generate sample of given `size` from geometric object
according to a homogeneous density. Optionally, provide `weights`
to specify custom sampling weights for the elements of a domain.
"""
struct LowDiscrepancySampling{T} <: ContinuousSamplingMethod
  nᵤ::Int
  dₘ::Int

  function LowDiscrepancySampling(nᵤ::Int, dₘ::Float64) where {T<:Real}
    if nᵤ ≤ 0
      throw(ArgumentError("Size must be positive"))
    end
    new{T}(nᵤ, dₘ)
  end
end

LowDiscrepancySampling(nᵤ::Int) = LowDiscrepancySampling(nᵤ, 1/√nᵤ)

LowDiscrepancySampling(dₘ::Float64) = LowDiscrepancySampling(round(1/dₘ²), dₘ)

function Sample(rng::AbstractRNG, mesh::SimpleMesh, method::LowDiscrepancySampling)
  for (idx, geom) in enumerate(mesh)
    if !isconvex(geom)
      throw(ArgumentError("Low discrepancy sampling only defined for convex geometries. 
      Item (#id: $idx) in the mesh is not a convex polytope."))
    end
  end
end

function rotate_shape(element::Triangle, floor::Bool = false)
  z_axis = Vec(0.0, 0.0, 1.0)  # Target normal (z-axis)
  rotation = Rotate(normal(element), z_axis)
  element  = element |> rotation
  if floor
      element = Triangle([Point(coords(point).x, coords(point).y, coords(point).z - coords(point).z) for point in vertices(element)]...)
  end
  element
end


function check_geometries(mesh::SimpleMesh)
  for element in mesh
      if !isa(element, Triangle)
          throw(ErrorException("This is not a triangle $element"))
      end
  end
end


function distance(p::Point{𝔼{3}}, q::Point{𝔼{3}})
  sqrt(
      (coords(p).x - coords(q).x)^2 +
      (coords(p).y - coords(q).y)^2 +
      (coords(p).z - coords(q).z)^2
  )
end


function distance(p::Point{𝔼{2}}, q::Point{𝔼{2}})
  sqrt(
      (coords(p).x - coords(q).x)^2 +
      (coords(p).y - coords(q).y)^2
  )
end


function get_longest_side(element::Triangle)
  tres = [(1,2),(2,3),(3,1)]
  vertices = vertices(element)
  longest_idx = 1

  for (idx, pair) in enumerate(tres)
      p1 = vertices[pair[1]]
      p2 = vertices[pair[2]]
      if distance(vertices[tres[longest_idx][1]], vertices[tres[longest_idx][1]]) < distance(p1, p2)
          longest_idx = idx
      end
  end
  tres[longest_idx]
end




# using LinearAlgebra
# using GeometryBasics

# # Function to compute the normal vector of a triangle
# # Input: `vertices` is a 3x3 matrix where each column is a vertex's [x, y, z] coordinates
# # Output: A unit vector normal to the plane of the triangle
# define_triangle_normal(vertices) = normalize(cross(
#     vertices[:, 2] - vertices[:, 1], # Vector from the first vertex to the second vertex
#     vertices[:, 3] - vertices[:, 1]  # Vector from the first vertex to the third vertex
# ))

# # Function to unfold a mesh to 2D
# # Inputs:
# #   `vertices`: A matrix where each column represents the [x, y, z] coordinates of a vertex in 3D
# #   `triangles`: A vector of tuples, where each tuple contains three indices defining a triangle
# # Output: A dictionary mapping vertex indices to their 2D positions after unfolding
# function unfold_mesh(vertices::Matrix{Float64}, triangles::Vector{NTuple{3, Int}})
#     # Map to store unfolded positions of vertices
#     unfolded_vertices = Dict{Int, Vector{Float64}}(undef, size(vertices, 2))

#     # Ensure the triangles list is not empty
#     if isempty(triangles)
#         error("The triangles input must not be empty.")
#     end

#     # Select the first triangle as the base, set its vertices in the 2D plane
#     base_triangle = triangles[1] # The first triangle in the mesh

#     # Helper function to validate a triangle index
#     function validate_index(idx::Int, max_idx::Int)
#         if idx < 1 || idx > max_idx
#             error("Triangle index $idx is out of bounds for the vertices array.")
#         end
#     end

#     # Validate triangle indices are within bounds of vertices
#     for idx in base_triangle
#         validate_index(idx, size(vertices, 2))
#     end

#     # Place the first vertex of the base triangle at the origin in 2D
#     unfolded_vertices[base_triangle[1]] = [0.0, 0.0]

#     # Place the second vertex along the x-axis at a distance equal to the edge length
#     unfolded_vertices[base_triangle[2]] = [
#         norm(vertices[:, base_triangle[2]] - vertices[:, base_triangle[1]]), # Edge length between vertex 1 and 2
#         0.0
#     ]

#     return unfolded_vertices
# end