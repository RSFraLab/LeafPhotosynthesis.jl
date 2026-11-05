"""
    LeafModel(; photosynthesis)

Small leaf composition boundary. The first implementation has one component,
`photosynthesis`, stored with its concrete type. Future components must stay
within leaf physiology; canopy radiation, air transport, soils, time stepping,
and device allocation belong outside this package.
"""
Base.@kwdef struct LeafModel{P<:AbstractPhotosynthesis}
    photosynthesis::P
end

@inline photosynthesize(leaf::LeafModel, env::LeafEnvironment) =
    photosynthesize(leaf.photosynthesis, env)
