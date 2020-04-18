function example_reorient(rootpath)
% example_reorient(rootpath)
% Simple example script that reorients all files found in a folder
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along
%    with this program; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    cd(rootpath)
    fileslist = ls(rootpath)
    fileslist = fileslist(3:end, :)  % remove '.' and '..' from the list
    auto_acpc_reorient(fileslist)
end %endfunction