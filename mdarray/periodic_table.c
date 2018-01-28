/***************************************************************************

    mdarray

    Python module for manipulation of atomic coordinates
    Copyright (C) 2012, Borys Szefczyk

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************/

#include "periodic_table.h"

/* Covalent radii are taken from:
 * http://pubs.rsc.org/en/content/articlehtml/2008/dt/b801115j
 */

Element element_table[] = {
{   1,   1.008      ,   "H",      "Hydrogen",  0.31 },
{   2,   4.002602   ,  "He",        "Helium",  0.28 },
{   3,   6.94       ,  "Li",       "Lithium",  1.28 },
{   4,   9.012182   ,  "Be",     "Beryllium",  0.96 },
{   5,  10.81       ,   "B",         "Boron",  0.84 },
{   6,  12.011      ,   "C",        "Carbon",  0.73 },
{   7,  14.007      ,   "N",      "Nitrogen",  0.71 },
{   8,  15.999      ,   "O",        "Oxygen",  0.66 },
{   9,  18.9984032  ,   "F",      "Fluorine",  0.57 },
{  10,  20.1797     ,  "Ne",          "Neon",  0.58 },
{  11,  22.98976928 ,  "Na",        "Sodium",  1.66 },
{  12,  24.3050     ,  "Mg",     "Magnesium",  1.41 },
{  13,  26.9815386  ,  "Al",     "Aluminium",  1.21 },
{  14,  28.085      ,  "Si",       "Silicon",  1.11 },
{  15,  30.973762   ,   "P",    "Phosphorus",  1.07 },
{  16,  32.06       ,   "S",        "Sulfur",  1.05 },
{  17,  35.45       ,  "Cl",      "Chlorine",  1.02 },
{  18,  39.948      ,  "Ar",         "Argon",  1.06 },
{  19,  39.0983     ,   "K",     "Potassium",  2.03 },
{  20,  40.078      ,  "Ca",       "Calcium",  1.76 },
{  21,  44.955912   ,  "Sc",      "Scandium",  1.70 },
{  22,  47.867      ,  "Ti",      "Titanium",  1.60 },
{  23,  50.9415     ,   "V",      "Vanadium",  1.53 },
{  24,  51.9961     ,  "Cr",      "Chromium",  1.39 },
{  25,  54.938045   ,  "Mn",     "Manganese",  1.50 },
{  26,  55.845      ,  "Fe",          "Iron",  1.42 },
{  27,  58.933195   ,  "Co",        "Cobalt",  1.38 },
{  28,  58.6934     ,  "Ni",        "Nickel",  1.24 },
{  29,  63.546      ,  "Cu",        "Copper",  1.32 },
{  30,  65.38       ,  "Zn",          "Zinc",  1.22 },
{  31,  69.723      ,  "Ga",       "Gallium",  1.22 },
{  32,  72.63       ,  "Ge",     "Germanium",  1.20 },
{  33,  74.92160    ,  "As",       "Arsenic",  1.19 },
{  34,  78.96       ,  "Se",      "Selenium",  1.20 },
{  35,  79.904      ,  "Br",       "Bromine",  1.20 },
{  36,  83.798      ,  "Kr",       "Krypton",  1.16 },
{  37,  85.4678     ,  "Rb",      "Rubidium",  2.20 },
{  38,  87.62       ,  "Sr",     "Strontium",  1.95 },
{  39,  88.90585    ,   "Y",       "Yttrium",  1.90 },
{  40,  91.224      ,  "Zr",     "Zirconium",  1.75 },
{  41,  92.90638    ,  "Nb",       "Niobium",  1.64 },
{  42,  95.96       ,  "Mo",    "Molybdenum",  1.54 },
{  43,  -1.0        ,  "Tc",    "Technetium",  1.47 },
{  44, 101.07       ,  "Ru",     "Ruthenium",  1.46 },
{  45, 102.90550    ,  "Rh",       "Rhodium",  1.42 },
{  46, 106.42       ,  "Pd",     "Palladium",  1.39 },
{  47, 107.8682     ,  "Ag",        "Silver",  1.45 },
{  48, 112.411      ,  "Cd",       "Cadmium",  1.44 },
{  49, 114.818      ,  "In",        "Indium",  1.42 },
{  50, 118.710      ,  "Sn",           "Tin",  1.39 },
{  51, 121.760      ,  "Sb",      "Antimony",  1.39 },
{  52, 127.60       ,  "Te",     "Tellurium",  1.38 },
{  53, 126.90447    ,   "I",        "Iodine",  1.39 },
{  54, 131.293      ,  "Xe",         "Xenon",  1.40 },
{  55, 132.9054519  ,  "Cs",       "Caesium",  2.44 },
{  56, 137.327      ,  "Ba",        "Barium",  2.15 },
{  57, 138.90547    ,  "La",     "Lanthanum",  2.07 },
{  58, 140.116      ,  "Ce",        "Cerium",  2.04 },
{  59, 140.90765    ,  "Pr",  "Praseodymium",  2.03 },
{  60, 144.242      ,  "Nd",     "Neodymium",  2.01 },
{  61,  -1.0        ,  "Pm",    "Promethium",  1.99 },
{  62, 150.36       ,  "Sm",      "Samarium",  1.98 },
{  63, 151.964      ,  "Eu",      "Europium",  1.98 },
{  64, 157.25       ,  "Gd",    "Gadolinium",  1.96 },
{  65, 158.92535    ,  "Tb",       "Terbium",  1.94 },
{  66, 162.500      ,  "Dy",    "Dysprosium",  1.92 },
{  67, 164.93032    ,  "Ho",       "Holmium",  1.92 },
{  68, 167.259      ,  "Er",        "Erbium",  1.89 },
{  69, 168.93421    ,  "Tm",       "Thulium",  1.90 },
{  70, 173.054      ,  "Yb",     "Ytterbium",  1.87 },
{  71, 174.9668     ,  "Lu",      "Lutetium",  1.87 },
{  72, 178.49       ,  "Hf",       "Hafnium",  1.75 },
{  73, 180.94788    ,  "Ta",      "Tantalum",  1.70 },
{  74, 183.84       ,   "W",      "Tungsten",  1.62 },
{  75, 186.207      ,  "Re",       "Rhenium",  1.51 },
{  76, 190.23       ,  "Os",        "Osmium",  1.44 },
{  77, 192.217      ,  "Ir",       "Iridium",  1.41 },
{  78, 195.084      ,  "Pt",      "Platinum",  1.36 },
{  79, 196.966569   ,  "Au",          "Gold",  1.36 },
{  80, 200.59       ,  "Hg",       "Mercury",  1.32 },
{  81, 204.38       ,  "Tl",      "Thallium",  1.45 },
{  82, 207.2        ,  "Pb",          "Lead",  1.46 },
{  83, 208.98040    ,  "Bi",       "Bismuth",  1.48 },
{  84,  -1.0        ,  "Po",      "Polonium",  1.40 },
{  85,  -1.0        ,  "At",      "Astatine",  1.50 },
{  86,  -1.0        ,  "Rn",         "Radon",  1.50 },
{  87,  -1.0        ,  "Fr",      "Francium",  2.60 },
{  88,  -1.0        ,  "Ra",        "Radium",  2.21 },
{  89,  -1.0        ,  "Ac",      "Actinium",  2.15 },
{  90, 232.03806    ,  "Th",       "Thorium",  2.06 },
{  91, 231.03588    ,  "Pa",  "Protactinium",  2.00 },
{  92, 238.02891    ,   "U",       "Uranium",  1.96 },
{  93,  -1.0        ,  "Np",     "Neptunium",  1.90 },
{  94,  -1.0        ,  "Pu",     "Plutonium",  1.87 },
{  95,  -1.0        ,  "Am",     "Americium",  1.80 },
{  96,  -1.0        ,  "Cm",        "Curium",  1.69 },
{  97,  -1.0        ,  "Bk",     "Berkelium", -1.0 },
{  98,  -1.0        ,  "Cf",   "Californium", -1.0 },
{  99,  -1.0        ,  "Es",   "Einsteinium", -1.0 },
{ 100,  -1.0        ,  "Fm",       "Fermium", -1.0 },
{ 101,  -1.0        ,  "Md",   "Mendelevium", -1.0 },
{ 102,  -1.0        ,  "No",      "Nobelium", -1.0 },
{ 103,  -1.0        ,  "Lr",    "Lawrencium", -1.0 },
{ 104,  -1.0        ,  "Rf", "Rutherfordium", -1.0 },
{ 105,  -1.0        ,  "Db",       "Dubnium", -1.0 },
{ 106,  -1.0        ,  "Sg",    "Seaborgium", -1.0 },
{ 107,  -1.0        ,  "Bh",       "Bohrium", -1.0 },
{ 108,  -1.0        ,  "Hs",       "Hassium", -1.0 },
{ 109,  -1.0        ,  "Mt",    "Meitnerium", -1.0 },
{ 110,  -1.0        ,  "Ds",  "Darmstadtium", -1.0 },
{ 111,  -1.0        ,  "Rg",   "Roentgenium", -1.0 },
{ 112,  -1.0        ,  "Cn",   "Copernicium", -1.0 },
{ 113,  -1.0        ,  "Nh",      "Nihonium", -1.0 },
{ 114,  -1.0        ,  "Fl",     "Flerovium", -1.0 },
{ 115,  -1.0        ,  "Mc",     "Moscovium", -1.0 },
{ 116,  -1.0        ,  "Lv",   "Livermorium", -1.0 },
{ 117,  -1.0        ,  "Ts",    "Tennessine", -1.0 },
{ 118,  -1.0        ,  "Og",     "Oganesson", -1.0 },
/*           SENTINEL           */
{  -1,  -1.0        ,  NULL,              "", -1.0 },
};


