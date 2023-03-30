#!/usr/bin/env python3
"""
Copyright 2020-2023 Edgardo M. Ortiz (e.ortiz.v@gmail.com)
https://github.com/edgardomortiz/Captus

This file is a convenience wrapper for running Captus assemble directly from the source tree. By
executing `captus_assemble-runner.py`, users can run captus_assemble without installing Captus.

This file is part of Captus. Captus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Captus is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Captus. If
not, see <http://www.gnu.org/licenses/>.
"""


from captus.captus_assembly import main


if __name__ == '__main__':
    main()
