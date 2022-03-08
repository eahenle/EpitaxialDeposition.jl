### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 185fce88-b4c3-4be7-84b4-7f2ce3db1c06
using DataFrames

# ╔═╡ bdf2b950-8afe-11ec-1732-61696ed36c9b
md"
## CHE 540 - Chemical Reactions
### Final Project Notebook

"

# ╔═╡ 82f306b3-0039-406e-a99a-68303ecdcfac
md"
## Source 1: [\"Lecture 08, CVD\"](https://research.ece.ncsu.edu/nanoengineering/files/Lecture%2008,%20CVD.ppt)
"

# ╔═╡ b256ab46-4990-4e6a-8966-453790fa774b
md"
### **Reactors and Systems** 
There are two main ways of producing thin films:
- Physical vapor deposition (PVD)
- Chemical vapor deposition (CVD)

In CVD, volatile, gaseous precursors are the source of the film species. There are a few common reactor configurations used in epitaxial Si deposition, which are \"horizontal\", \"barrel\", and \"pancake\".

*insert image*

There are also many types of CVD, such as 
- *Ambient pressure* (AP) CVD 
- *Low pressure* (LP) CVD - typically ranges from **10 mT to 1 T**, to increase diffusion and the mean free path of gases in moving toward the substrate
- *Plasma enhanced* (PE) CVD - lowers overall reaction temperature, as the precursors are heated to plasma, which react quickly with the substrate surface
- *Vapor phase epitaxy* (VPE) - having a single crystal subtrate act as a template for a film of identical or similar crystal structure

### **Chemistry**
Both heterogeneous and homogeneous reactions are part of the CVD process.

*insert image*

These steps include the fundamental steps of surface catalyzed reactions which have been discussed in class before:
1) Transport of reactant species to the boundary layer of the solid
2) Chemical reaction between the reactants
3) Adsorption onto the material's surface
4) Surface kinetics, growth of the film
5) Desorption from the surface
6) Transport from the boundary layer to the bulk source, out of the process

Wafers used in the semiconductor industry are produced via **VPE** where there are multiple possible reactants, such as silicon tetrachloride (``SiCl_4``) or silicon hydride (a silane - ``SiH_4``).

- 1200 °C: `` SiCl_4 (g) + 2H_2 (g) \rightarrow Si(s) + 4HCl (g)`` 
- 650 °C: ``SiH_4 (g) \rightarrow Si(s) + 2H_2 (g)``


*insert image* 


"

# ╔═╡ 0ae4bec3-f8dd-4c99-b168-a3d48b113dda
md"
### Source 2: [\"Modelling of Si VPE using SiCl₄ as a Source \"](https://www.sciencedirect.com/science/article/pii/0026269295000127)
"

# ╔═╡ 60fccd3e-e87b-49b5-b073-dc16dc4e0c92
md"
## **Mathematical Modeling**


This paper proposes a growth-rate model based on chemical kinetics for VPE of silicon, using **SiCl₄** as a source, in a **horizontal rectangular reactor at atmospheric pressure**. The model includes factors such as **temperature, flow rate, mole fraction, and position,** and the model agrees well with experimental data available.

A few resources are referenced in this paper which seem to cover epitaxy at low pressures:
- \[11\] *K.F. Jensen and D.B. Graves, Modelling and analysis of low pressure CVD reactor, J. Electrochem. Soc., 131 (1983) 1950.*
- \[12\] *M.G. Joshi, Modelling of LPCVD reactors, J. Electrochem. Soc., 134 (1987) 318.*

A simplified reaction scheme is described:
1. ``SiCl_4 + 2H_2 \rightarrow Si (s) + 4HCl``, deposition of Si on the substrate surface
2. ``Si + 2HCl \rightarrow SiCl_2 + H_2``, etching of Si by HCl
3. ``SiCl_4 + Si (s) \rightarrow 2SiCl_2``, an additional etching reaction, when the concentration of SiCl₄ is in the mixture is large

<s>We note an additional reaction set which deals with the extremely unstable and short-lived SiCl₂:
4. ``SiCl_2 + H_2 \rightarrow SiCl_2H_2``<\s>

"

# ╔═╡ e2678619-a552-4397-af1b-f657c03c58a8
md"
## **Article Nomenclature and Kinetic Data**

The final models for total silicon deposition rate, *F*, and total silicon film growth rate, *G*, are given in the paper.

1. **Silicon deposition rate:** ``F = (H(C_{BS}-C_{SS})-K_3C_{SS}) \left[1-\frac{2K_2}{H_{HCl}+K_2} \right] - K_3C_{SS}``

* ``H`` = gas phase mass transfer coefficient of SiCl₄, ``D(\frac{\eta X}{\rho V})^{-1/2}``
  * D = diffusivity of SiCl₄ in H₂
  * ``\eta`` = viscosity of the carrier gas
  * ``\rho`` = density of the carrier gas
  * ``X`` = the point along the x-direction in the (rectangular) reactor

* ``H_{HCl}`` = gas phase mass transfer coefficient of HCl, the value of which is reported to be ~3x that of SiCl₄
* ``C_{BS}`` = bulk concentration of SiCl₄
* ``C_{SS}`` = concentration of SiCl₄ in the gas phase adjacent to growing film


* ``K`` = the rate constant for reaction (1), the surface reaction of Si with H₂, ``K = K_0exp[\frac{-E}{RT}]``
* ``K_2`` = the rate constant for reaction (2), the surface reaction of HCl with Si, ``K_2 = K_{20}exp[\frac{-E_{a2}}{RT}]``
* ``K_3`` = the rate constant for reaction (3), the surface reaction of SiCl₄ on Si, ``K_3 = K_{30}exp[\frac{-E_{a3}}{RT}]``

2. **Silicon film growth rate:** ``G = \frac{F}{N}``

* ``N`` = number of silicon atoms incorporated in a unit volume of the film

Values of the various reaction constants presented are not available in literature were estimated numerically by the authors to be:

* ``K_{0} = 2950.94 \frac{cm}{s}``, ``E = 18.910 \frac{kcal}{mol}``
* ``K_{20} = 39.3886 \frac{cm}{s}``, ``E_{a2} = 17.490 \frac{kcal}{mol}``
* ``K_{30} = 952.058 \frac{cm}{s}``, ``E_{a3} = 21.721 \frac{kcal}{mol}``

"

# ╔═╡ f58b39cc-3a7c-4b30-bda5-6c837c2450cd
begin
	article_symbols = [:F, :G]
	
	definitions = [
		md"Total rate of silicon deposition on the substrate", 
		md"Total rate of silicon film growth on the substrate", 
	]
		# md"Bulk concentration of SiCl₄", 
		# md"Concentration of SiCl₄ in the gas phase adjacent to the growing surface"
	math_definitions = [
		md"$F_{1} = (H(C_{BS}-C_{SS}) - K_3C_{SS}) \left[ 1- \frac{2K_2}{H_{HCl}+K_2} \right]$",
		md"$D{{\left( \frac{\eta X} {\rho V} \right)}^{-1/2}}$"
	]
	
	DataFrame("Article Symbols" => article_symbols, 
			"Definition" => definitions, 
			"Symbolics" => math_definitions)
	
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"

[compat]
DataFrames = "~1.3.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═185fce88-b4c3-4be7-84b4-7f2ce3db1c06
# ╟─bdf2b950-8afe-11ec-1732-61696ed36c9b
# ╟─82f306b3-0039-406e-a99a-68303ecdcfac
# ╟─b256ab46-4990-4e6a-8966-453790fa774b
# ╟─0ae4bec3-f8dd-4c99-b168-a3d48b113dda
# ╟─60fccd3e-e87b-49b5-b073-dc16dc4e0c92
# ╠═e2678619-a552-4397-af1b-f657c03c58a8
# ╠═f58b39cc-3a7c-4b30-bda5-6c837c2450cd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
