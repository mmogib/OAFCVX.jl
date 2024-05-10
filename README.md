# OAFCVX

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmogib.github.io/OAFCVX.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmogib.github.io/OAFCVX.jl/dev/)
[![Build Status](https://github.com/mmogib/OAFCVX.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mmogib/OAFCVX.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Build Status](https://app.travis-ci.com/mmogib/OAFCVX.jl.svg?branch=master)](https://app.travis-ci.com/mmogib/OAFCVX.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mmogib/OAFCVX.jl?svg=true)](https://ci.appveyor.com/project/mmogib/OAFCVX-jl)
[![Coverage](https://codecov.io/gh/mmogib/OAFCVX.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmogib/OAFCVX.jl)

# Notes

## Directions

1. Fixed: d
2. Fixed p: p in P, d = p-v
3. Ideal Point: IP. solve (v-yI +epsilon)^T d =1
4. Adj (not doing it)
5. AAdj (new) based on how close to v and follow Adj procedure

## Vertices

1. Vertex selection with clusters (C)
2. Vertex selection with adjacency information (Adj)
3. Approximate Vertex selection with adjacency information (AAdj)
4. Vertex selection using local upper bounds (UB)
