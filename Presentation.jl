### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

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

**Volatile, gaseous precursors** are the source of the species to be deposited. There are a few common reactor configurations used in epitaxial Si deposition, but here we consider a rectangular reactor.

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

Vapor phase epitaxy is one technique for generating highly isotropic films, and there are multiple reactant feedstocks used, such as **silicon tetrachloride** (``SiCl_4``) or silicon hydride (a silane - ``SiH_4``).

- 1200 °C: `` SiCl_4 (g) + 2H_2 (g) \rightarrow Si(s) + 4HCl (g)`` 
- 650 °C: ``SiH_4 (g) \rightarrow Si(s) + 2H_2 (g)``

These reactions have a few nasty byproducts, including **hydrogen chloride** and **silicon dichloride**, which is a highly unstable molecule.
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

These are composite rate constants which include both flux to the surface, as well as the surface reaction itself, and their values:

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
4. **Lennard-Jones parameters** (molecular collision diameter) are estimated from each species' **critical temperatures and pressures**
    - SiCl₂ critical properties are not available and a quick DFT calculation attempting to elucidate them failed, so we used SiH₂Cl₂ as a proxy molecule
5. As outlined in the paper we based the model off, these **kinetic parameters are most valid around 1200 °C**
6. **Dilute conditions** for our limiting reactant and excess of others; **ideal gas behavior** of these species
7. Gas film at the surface of the wafer is **homogeneous**, allowing us to run a **single kinetic simulation** across the surface of the wafer 
8. Reactions only happen on the silicon surface

"""

# ╔═╡ 7a941657-7ae6-4eaf-8e06-7a6271a47a2d
md"""
## **Numerical Methods and Boundary/Initial Conditions**

⭐ This project has been done entirely in Julia -- something which we didn't expect to be able to say at the beginning. 

1. The general form of the partial differential equations that we were able to solve in Julia using a finite volume method is as follows: 


$$\alpha \frac{\partial \phi}{\partial t} + \nabla \cdot(-\mathscr{D} \nabla \phi) + \beta \phi = \gamma$$ 

- In the case of diffusion, this equation becomes:

$$\frac{\partial c}{\partial t} + \nabla \cdot(-\mathscr{D} \nabla c) = 0$$

2. The corresponding boundary conditions for this PDE (both in simplified and unsimplified form) are:

$$a \nabla \phi \cdot \textbf{n} + b \phi = c$$

- In this simulation we used **Neumann boundary conditions**, allowing us to specify species' flux at any cell boundary. This made it possible to create a symbolic "wafer" which spans the finite volume cells in the diffusion portion of the simulation

3. Initial conditions relating to species concentrations are set before the simulation, and are interactive / easily manipulated


4. Symbolic species ``Si_{dep}`` and ``Si_{etch}`` are generated as "products" in our reactions, and the difference between them is the total amount of silicon deposited across the wafer 

- The differential equations for the species are all in the gas phase, and were reaching unrealistic equilibria (not seeing depletion of reactants over time)


5. Film growth rate over time was calculated--since we had such small timesteps--as the change in film thickness over the timestep

"""

# ╔═╡ 695d2cd2-41dc-4749-a89e-fa187667bf47
md"""
## **Julia Packages and Acknowledgements**

- *Catalyst.jl* - allowed for the painless generation of the chemical reaction network and associated differential equations
- *JFVM.jl* - provided a framework for finite volume methods, allowing us to include diffusion of species from the wafer surface
- *PlutoUI.jl* - provided numerous tools for enabling interactivity within the notebook
- *LsqFit.jl* - used to fit a model to tabulated data from Welty, pertinent to the diffusion coefficients

🍮 **Dr. Kelsey Stoerzinger**, for helpful discussion and project guidance 🍮
"""

# ╔═╡ 840fa195-ee05-41a7-a2f3-6ac04b538ccb
md"""
## **Future Work**

- Realistic reactant sources and an outlet for the byproduct gas
- Finding **relevant experimental data**, which refines or extends the current model
- **Simultaneous reaction and diffusion** vs. two solvers which communicate between one another (+ fluid flow if possible)
- Introduction of **dopants** to create **P/N-type semiconductors**
- Different **reactor setups/geometries**, **multi-wafer growth**
- Finding more advanced reaction mechanisms

"""

# ╔═╡ a50df40f-0a00-456d-ba70-b5ecc637e6b4
md"""
## **Other References & Resources**

- **Project inspiration:** *Modelling of silicon epitaxy using silicon tetrachloride as the source*, D.K. Pal, M.K. Kowar, A.N. Daw, P. Roy, 1995

- **Spherical molecule diffusion model:** *Fundamentals of Momentum, Heat, and Mass Transfer, 6th Ed*, James R. Welty, Gregoryb L. Rorrer, David G. Foster, 2015

- **SiH₂Cl₂ critical temperature and pressure:** *Process Feasibility Study in Support of Silicon Materials*, Ku-Yen Li, Keith C. Hansen, Carl L. Yaws, June 1979

- **Equilibrium rate constant data for ``SiCl_4 + Si(s) \leftrightarrow SiCl_2``:** *Experimentelle Untersuchung des Reaktionsgleichgewichtes SiCl₄(g) + Si(f) = 2SiCl₂(g) nach der Stromungshmethode*, Von R. Teichmann, E. Wolf
"""

# ╔═╡ Cell order:
# ╟─bdf2b950-8afe-11ec-1732-61696ed36c9b
# ╟─82f306b3-0039-406e-a99a-68303ecdcfac
# ╟─b256ab46-4990-4e6a-8966-453790fa774b
# ╟─0ae4bec3-f8dd-4c99-b168-a3d48b113dda
# ╟─60fccd3e-e87b-49b5-b073-dc16dc4e0c92
# ╟─e2678619-a552-4397-af1b-f657c03c58a8
# ╟─b29aeadf-4f5a-42a3-942e-22559a7cab21
# ╠═7a941657-7ae6-4eaf-8e06-7a6271a47a2d
# ╠═695d2cd2-41dc-4749-a89e-fa187667bf47
# ╠═840fa195-ee05-41a7-a2f3-6ac04b538ccb
# ╠═a50df40f-0a00-456d-ba70-b5ecc637e6b4
