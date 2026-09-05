# Writes a `CGECameroon.RawData` out in `load_data`'s generic workbook layout (see
# DATA.md §5 and `load_data`'s docstring): sheets `sectors`, `labour`, `iotable`, `imat`,
# `employment`, `wagedist`, `miscellaneous`, `scalars`.
#
# This is the *test-side* reference implementation of that format -- the production writer
# is an R exporter -- and it exists so `runtests.jl` can (a) round-trip the real Cameroon
# workbook through the generic path and (b) build a synthetic N-sector dataset without
# checking a second binary .xlsx into the repo.
#
# Used by test/runtests.jl; not loaded by cge.jl.

using XLSX

"Writes one labelled matrix sheet: A1 a free corner, row 1 the column headers, column A
the row labels, `d[rowkey, colkey]` in the body."
function write_matrix_sheet!(xf, name::AbstractString,
                             rowkeys::Vector{Symbol}, colkeys::Vector{Symbol},
                             d::Dict{Tuple{Symbol,Symbol},Float64})
    s = XLSX.addsheet!(xf, name)
    s[1, 1] = name                       # free corner -- the loader ignores it
    for (j, c) in enumerate(colkeys)
        s[1, j + 1] = String(c)
    end
    for (i, r) in enumerate(rowkeys)
        s[i + 1, 1] = String(r)
        for (j, c) in enumerate(colkeys)
            s[i + 1, j + 1] = d[r, c]
        end
    end
    return s
end

"""
    write_generic_workbook(path, raw::RawData; sector_order, labour_order, labels...)

Writes `raw` to `path` in the generic layout, so that `load_data(path)` returns the same
sets and matrices. `sector_order`/`labour_order` control the order the codes appear in the
file (the loader reads by header lookup, so any permutation must load identically) --
`raw.IT` still decides each sector's `traded` flag.
"""
function write_generic_workbook(path::AbstractString, raw;
                                sector_order::Vector{Symbol} = raw.SEC,
                                labour_order::Vector{Symbol} = raw.LC,
                                sector_labels::AbstractDict = Dict{Symbol,String}(),
                                labour_labels::AbstractDict = Dict{Symbol,String}())
    isfile(path) && rm(path)
    XLSX.openxlsx(path, mode = "w") do xf
        sec = xf[1]
        XLSX.rename!(sec, "sectors")
        sec[1, 1] = "code"; sec[1, 2] = "label"; sec[1, 3] = "traded"
        for (n, i) in enumerate(sector_order)
            sec[n + 1, 1] = String(i)
            sec[n + 1, 2] = get(sector_labels, i, String(i))
            sec[n + 1, 3] = i in raw.IT
        end

        lab = XLSX.addsheet!(xf, "labour")
        lab[1, 1] = "code"; lab[1, 2] = "label"
        for (n, l) in enumerate(labour_order)
            lab[n + 1, 1] = String(l)
            lab[n + 1, 2] = get(labour_labels, l, String(l))
        end

        write_matrix_sheet!(xf, "iotable",       sector_order, sector_order, raw.io)
        write_matrix_sheet!(xf, "imat",          sector_order, sector_order, raw.imat)
        write_matrix_sheet!(xf, "employment",    sector_order, labour_order, raw.xle)
        write_matrix_sheet!(xf, "wagedist",      sector_order, labour_order, raw.wdist)
        write_matrix_sheet!(xf, "miscellaneous", CGECameroon.MISC_ROWS, sector_order, raw.zz)

        sca = XLSX.addsheet!(xf, "scalars")
        sca[1, 1] = "name"; sca[1, 2] = "value"
        r = 1
        for name in CGECameroon.SCALAR_NAMES
            haskey(raw.scalars, name) || continue
            r += 1
            sca[r, 1] = String(name)
            sca[r, 2] = raw.scalars[name]
        end
        for l in labour_order
            r += 1
            sca[r, 1] = "wa0_" * String(l)
            sca[r, 2] = raw.wa0[l]
        end
    end
    return path
end
