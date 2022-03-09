### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 185fce88-b4c3-4be7-84b4-7f2ce3db1c06
using DataFrames, PlutoUI

# ╔═╡ bdf2b950-8afe-11ec-1732-61696ed36c9b
md"
# CHE 540 - Chemical Reactions I 
## **Final Project Presentation Notebook**
### Adrian Henle, Ian Harreschou

\
⭐ Welcome everyone, and thank you for coming to our presentation! Today we'll be walking you through our attempt at modeling a **chemical vapor deposition** (CVD) reactor, more specifically one used for **silicon vapor phase epitaxy** (Si VPE).

We hope to give you a brief yet in-depth tour of our simplifying assumptions, methods, and most excitingly the current reactor model!
"

# ╔═╡ 82f306b3-0039-406e-a99a-68303ecdcfac
md"
### Source 1: [\"Lecture 08, CVD\"](https://research.ece.ncsu.edu/nanoengineering/files/Lecture%2008,%20CVD.ppt)
"

# ╔═╡ b256ab46-4990-4e6a-8966-453790fa774b
md"
### **Reactors and Systems** 
Chemical vapor deposition is a manufacturing technique used throughout the semiconductor industry to produce high performance thin films. 

In CVD, volatile, gaseous precursors are the source of the species to be deposited. There are a few common reactor configurations used in epitaxial Si deposition, wbut here we consider a rectangular reactor.

*insert image*

There are also many types of CVD, such as 
- *Ambient pressure* (APCVD) 
- *Low pressure* (LPCVD)- typically ranges from 10 mT to 1 T, to increase diffusion and the mean free path of gases in moving toward the substrate
- *Plasma enhanced* (PECVD) - lowers overall reaction temperature, as the precursors are heated to plasma, which react quickly with the substrate surface
- **Vapor phase epitaxy (VPE) - having a single crystal subtrate act as a template for a film of identical or similar crystal structure**

### **Chemistry**

Reaction steps follow closely to microkinetic models for surface reactions that we've seen in lectures to this point:
1) Transport of reactant species to the boundary layer of the solid
2) Chemical reaction between the reactants
3) Adsorption onto the material's surface
4) Surface kinetics, growth of the film
5) Desorption from the surface
6) Transport from the boundary layer to the bulk gas, and out of the reactor

Vapor phase epitaxy is one technique for generating highly isotropic films, and there are multiple reactant feedstocks used, such as silicon tetrachloride (``SiCl_4``) or silicon hydride (a silane - ``SiH_4``).

- 1200 °C: `` SiCl_4 (g) + 2H_2 (g) \rightarrow Si(s) + 4HCl (g)`` 
- 650 °C: ``SiH_4 (g) \rightarrow Si(s) + 2H_2 (g)``

These reactions have a few nasty byproducts, including **hydrogen chloride** and **silicon dichloride**, which is a highly unstable molecule.


*insert image* 


"

# ╔═╡ 0ae4bec3-f8dd-4c99-b168-a3d48b113dda
md"
### Source 2: [\"Modelling of Si VPE using SiCl₄ as a Source \"](https://www.sciencedirect.com/science/article/pii/0026269295000127)
"

# ╔═╡ 60fccd3e-e87b-49b5-b073-dc16dc4e0c92
md"
## **Mathematical Modeling**


⭐ Luckily, we were quickly able to find a paper dealing with exactly what we wanted to explore: numerical modeling of silicon epitaxy. 

This paper proposes a growth-rate model using **SiCl₄** as a source, in a **horizontal rectangular reactor at atmospheric pressure**. The model includes factors such as **temperature, flow rate, mole fraction, and position within the reactor,** and the model--after some numerical manipulation--agrees well with experimental data available.

!!! note
	A few resources are referenced in this paper which seem to cover epitaxy at low pressures:
	- \[11\] *K.F. Jensen and D.B. Graves, Modelling and analysis of low pressure CVD reactor, J. Electrochem. Soc., 131 (1983) 1950.*
	- \[12\] *M.G. Joshi, Modelling of LPCVD reactors, J. Electrochem. Soc., 134 (1987) 318.*

The authors list the following for a simplified reaction scheme:
1. ``SiCl_4 + 2H_2 \rightarrow Si (s) + 4HCl``, deposition of Si on the substrate surface
2. ``Si + 2HCl \rightarrow SiCl_2 + H_2``, etching of Si by HCl
3. ``SiCl_4 + Si (s) \rightarrow 2SiCl_2``, an additional etching reaction, when the concentration of SiCl₄ is in the mixture is large

The decomposition of SiCl₄ is known to occur **on the surface** *rather than in the gas phase*, leaving behind single crystal layers of silicon on the substrate, and this rate is proportional to the number of surface sites occupied by SiCl₄. A **Langmuir adsorption model** is assumed by this analysis, and is therefore baked into our assumptions.

"

# ╔═╡ e2678619-a552-4397-af1b-f657c03c58a8
md"
## **Article Nomenclature and Kinetic Data**


* ``K`` = the rate constant for reaction (1), the surface reaction of Si with H₂, ``K = K_0exp[\frac{-E}{RT}]``
* ``K_2`` = the rate constant for reaction (2), the surface reaction of HCl with Si, ``K_2 = K_{20}exp[\frac{-E_{a2}}{RT}]``
* ``K_3`` = the rate constant for reaction (3), the surface reaction of SiCl₄ on Si, ``K_3 = K_{30}exp[\frac{-E_{a3}}{RT}]``

These are composite rate constants which include both flux to the surface, as well as the surface reaction itself. Values of the various reaction constants presented were not available in literature (and relevant experimental data is extremely scarce it seems), so they were numerically fit by the authors to experimental data, coming out to be:

* ``K_{0} = 2950.94 \frac{cm}{s}``, ``E = 18.910 \frac{kcal}{mol}``
* ``K_{20} = 39.3886 \frac{cm}{s}``, ``E_{a2} = 17.490 \frac{kcal}{mol}``
* ``K_{30} = 952.058 \frac{cm}{s}``, ``E_{a3} = 21.721 \frac{kcal}{mol}``

⭐ In our analysis, these constants are multiplied by our wafer surface area to effectively scale the reaction rates (larger wafers will result in higher overall rate and quicker depletion of reactants) -- something that wasn't being considered prior!

"

# ╔═╡ b29aeadf-4f5a-42a3-942e-22559a7cab21
md"""
## **Simplifying Assumptions**

1. In the simplest model, **gaseous species are initially instantaneously equilibrated with the reactor volume**, resulting in no concentration gradient
2. The **inlet species are only SiCl₄ and H₂**, their mole fractions being complementary of one another
3. The **reactor operates isothermally**
4. **Lennard-Jones parameters** (molecular collision diameter) are estimated from each species' **critical temperatures and pressures** in calculation diffusion coefficients, the method for which is outlined in Ch. 24 of *Fundamentals of Momentum, Heat, and Mass Transfer, 6th Ed*
    - SiCl₂ critical properties are not available and a quick DFT calculation attempting to elucidate them failed, so we used SiH₂Cl₂ as a proxy molecule
5. As outlined in the paper we based the model off, these **kinetic parameters are most valid around 1200 °C**
6. **Dilute conditions** for our limiting reactant and excess of others; **ideal gas behavior** of these species
7. Gas film at the surface of the wafer is homogeneous, allowing us to run a single kinetic simulation across the surface of the wafer (opposed to independent simulations for each "cell" in our diffusion model)

"""

# ╔═╡ 7a941657-7ae6-4eaf-8e06-7a6271a47a2d
md"""
## **Numerical Methods and Boundary/Initial Conditions**

⭐ This project has been done entirely in Julia -- something which we didn't expect to be able to say at the beginning. 

The general form of the partial differential equations that we were able to solve in Julia using a finite volume method is as follows: 


$$\alpha \frac{\partial \phi}{\partial t} + \nabla \cdot(-\mathscr{D} \nabla \phi) + \beta \phi = \gamma$$

with $$\phi$$ being the variable of interest. In the case of diffusion, this equation becomes:

$$\frac{\partial c}{\partial t} + \nabla \cdot(-\mathscr{D} \nabla c) = 0$$

The corresponding boundary conditions for this PDE (both in simplified and unsimplified form) are:

$$a \nabla \phi \cdot \textbf{n} + b \phi = c$$

where in this simulation we used Neumann boundary conditions, allowing us to specify ∅ flux at any cell boundary.

initial conditions which relate to concentration of species (set during the simulation)

how we describe the "pseudo-concentration" of deposited silicon (Si_dep and Si_etch)

reactions only happen on the silicon surface

film growth rate over time, simple derivative of growth array over Δt

"""

# ╔═╡ 695d2cd2-41dc-4749-a89e-fa187667bf47
md"""
## **Julia Packages and Acknowledgements**


"""

# ╔═╡ 840fa195-ee05-41a7-a2f3-6ac04b538ccb
md"""
## **Future Work**

- Instead of assuming that our gases are instantly dosed into our reactor and equilibrate with its volume, we have an inlet nozzle for reactants and an outlet for the byproduct gas (via manipulating the boundary conditions)
- Accumulating more relevant experimental data, and hopefully data which corroborates the current model
- Complete model for simultaneous reaction and diffusion, rather than using two solvers which communicate between one another (+ fluid flow if possible)
- Introduction of dopants to create P/N-type semiconductors
- Investigation of different reactor setups, multi-wafer growth

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DataFrames = "~1.3.2"
PlutoUI = "~0.7.37"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

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

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

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
# ╟─e2678619-a552-4397-af1b-f657c03c58a8
# ╠═b29aeadf-4f5a-42a3-942e-22559a7cab21
# ╠═7a941657-7ae6-4eaf-8e06-7a6271a47a2d
# ╠═695d2cd2-41dc-4749-a89e-fa187667bf47
# ╠═840fa195-ee05-41a7-a2f3-6ac04b538ccb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
