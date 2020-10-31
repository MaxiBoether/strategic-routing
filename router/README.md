# Routing Framework
## Description
This C++ project offers a framework for research projects regarding (strategic)
routing algorithms. It is the result of our work on our ATMOS paper (doi:
10.4230/OASIcs.ATMOS.2020.10) as well as several bachelor theses at the chair
for Algorithm Engineering at the Hasso Plattner Insitute in Potsdam. For the
remainder of this readme, we assume familiarity with the ATMOS paper.

## Project Structure

The repository contains the outer framework as well as the implementation of our
algorithms solving SAP. You find all implementation sourcecode in the `src`
folder. The subdirectory `core` contains the core router part
(framework), responsible for input/output/general datatypes. The other
subdirectores are for each routing submodule ("strategy"). In this repository,
there only is the `ssotd` strategic, which maps to the SAP problem in the paper.
For example, we also implemented other evolutionary algorithms solving SAP,
which define a different strategy. Strategies can have different variants (e.g.
1D-SAP).

Used libraries are found in `lib`. When building, you can find the results in
`build`. All header files should be in `include`.

Please note that due to historic reasons, the naming of problems/algorithms in
the sourcecode differs to the naming in the paper. `SSOTD` maps to the `SAP`
problem in general. The `fulldisjoint` variant maps to `D-SAP`, the
`newnodisjoint` variant maps to `SAP`, and the `newonedisjoint` variant maps to
`1D-SAP`. If the `new` prefix on the variants is missing, that indicates the
fewer criteria variants (´FC`).

In the SAP implementations, you will find several optimization techniques not
yet described in any paper. We plan to publish the proofs of correctness in the
journal version of the paper.

## Input and Output Data Format
Due to historic reasons, we are using MATSim's graph and plan format. We refer
to the MATSim user guide for details (https://www.matsim.org/docs/userguide/).
Unfortunately, we cannot provide you with our graph representing on Berlin. The
idea in the `plans.xml` file is that you can specify multiple persons with the
same OD-pair (and starting time) to set how many people should be routed. SAP is
solved for each OD-pair. The framework outputs a plans.xml that can be executed
in the MATSim simulator.

This means you can validate your (and our) algorithms in simulation settings. We
recommend Simunto Via for visulization. Unfortunately, we did not have the time
in our bachelor project to work on simulation evaluations of our algorithm. This
would be a good continuation of our work.

## Building
### Requirements
Make sure to have the following libraries installed on your system.

- Tinyxml2
- GNU Scientific Library (GSL), including the BLAS
- OpenMP

### Makefile Usage
You can build the project using `make` in the project root. By default, we use some debug and protection flags:

```
-g -fstack-protector-strong -fstack-clash-protection -fcf-protection -Wall -Wpedantic -Wextra -std=c++2a -D_GLIBCXX_ASSERTIONS
```

Furthermore, you can specify `SANITIZE=thread` or `SANITIZE=address` to either include the thread or address sanitizer in your build (default is address). When specifying `TYPE=DEBUG`, in addition to the flags above, the `-O0` flag is added. Without `TYPE=DEBUG`, we use `-O2 -D_FORTIFY_SOURCE=2` instead.

When specifying TYPE=RELEASE, all of these flags are omitted and `-Ofast` is added instead.

To select the used psychological model, you can use the PSYCHMOD variable. We default to user_equilibrium_2r.

To select which module you are building, you can use the STRATEGY variable. We default to SSOTD fulldisjoint. Possible strategies are:
- sstod. In this case, you can also specify SSOTD_VARIANT, which can either be
  onedisjoint, nodisjoint, newonedisjoint, newnodisjoint or fulldisjoint, which is our default.
- The other strategies we worked on are not (yet) published, so we do not
  include them for now.

You can also add more strategies just by creating more subfolders in src. Please
remember to put your headers in include, as this is added to the include path.
Also make sure, if you need more cpp files than one, to expand the Makefile
accordingly. For this, you have to modify the ADDITIONALS variable. Just take a
look at the example for the ea or ssotd and I'm sure you can figure it out.

## Usage
```
./router <graph> <plans> <output>
```

### Libraries and Licenses
We employ some external libraries. We would like to thank the authors for their
work

Argh! A minimalist argument handler. (https://github.com/adishavit/argh)

Copyright (c) 2016, Adi Shavit
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of  nor the names of its contributors may be used to
   endorse or promote products derived from this software without specific
   prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


JSON for Modern C++ (https://github.com/nlohmann/json)

Licensed under the MIT License <http://opensource.org/licenses/MIT>.
SPDX-License-Identifier: MIT
Copyright (c) 2013-2019 Niels Lohmann <http://nlohmann.me>.
Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Tinyxml2 (https://github.com/leethomason/tinyxml2)

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any
damages arising from the use of this software.

Permission is granted to anyone to use this software for any
purpose, including commercial applications, and to alter it and
redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and
must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.

### Project License

This soure code is licensed under GNU General Public License v3 (GPL-3) (see
`LICENSE.md`).
