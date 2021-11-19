/*Copyright 2019 DScribe developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#include "celllist.h"
#include <algorithm>
#include <map>
#include <math.h>
#include <utility>

using namespace std;

CellList::CellList(py::array_t<double> positions, double cutoff)
    : positions(positions.unchecked<2>()), cutoff(cutoff),
      cutoffSquared(cutoff * cutoff) {
  if (cutoff > 0) {
    this->init();
  }
}

void CellList::init() {
  // Find cell limits
  this->xmin = this->xmax = this->positions(0, 0);
  this->ymin = this->ymax = this->positions(0, 1);
  this->zmin = this->zmax = this->positions(0, 2);
  for (ssize_t i = 0; i < this->positions.shape(0); i++) {
    double x = this->positions(i, 0);
    double y = this->positions(i, 1);
    double z = this->positions(i, 2);

    if (x < this->xmin) {
      this->xmin = x;
    } else if (x > this->xmax) {
      this->xmax = x;
    };

    if (y < this->ymin) {
      this->ymin = y;
    } else if (y > this->ymax) {
      this->ymax = y;
    };

    if (z < this->zmin) {
      this->zmin = z;
    } else if (z > this->zmax) {
      this->zmax = z;
    };
  };

  // Add small padding to avoid floating point precision problems at the
  // boundary
  constexpr double padding = 0.0001;
  this->xmin -= padding;
  this->xmax += padding;
  this->ymin -= padding;
  this->ymax += padding;
  this->zmin -= padding;
  this->zmax += padding;

  // Determine amount and size of bins. The bins are made to be always of equal
  // size.
  this->nx = max(1, int((this->xmax - this->xmin) / this->cutoff));
  this->ny = max(1, int((this->ymax - this->ymin) / this->cutoff));
  this->nz = max(1, int((this->zmax - this->zmin) / this->cutoff));
  this->dx = max(this->cutoff, (this->xmax - this->xmin) / this->nx);
  this->dy = max(this->cutoff, (this->ymax - this->ymin) / this->ny);
  this->dz = max(this->cutoff, (this->zmax - this->zmin) / this->nz);

  // Initialize the bin data structure. It is a 4D vector.
  this->bins = vector<vector<vector<vector<int>>>>(
      this->nx, vector<vector<vector<int>>>(
                    this->ny, vector<vector<int>>(this->nz, vector<int>())));

  // Fill the bins with atom indices
  for (ssize_t idx = 0; idx < this->positions.shape(0); idx++) {
    double x = this->positions(idx, 0);
    double y = this->positions(idx, 1);
    double z = this->positions(idx, 2);

    // Get bin index
    int i = (x - this->xmin) / this->dx;
    int j = (y - this->ymin) / this->dy;
    int k = (z - this->zmin) / this->dz;

    // Add atom index to the bin
    this->bins[i][j][k].push_back(idx);
  };
}

CellListResult CellList::getNeighboursForPosition(const double x,
                                                  const double y,
                                                  const double z) const {
  // The indices of the neighbouring atoms
  vector<int> neighbours;
  vector<double> distances;
  vector<double> distancesSquared;

  // Find bin for the given position
  int i0 = (x - this->xmin) / this->dx;
  int j0 = (y - this->ymin) / this->dy;
  int k0 = (z - this->zmin) / this->dz;

  // Get the bin ranges to check for each dimension.
  int istart = max(i0 - 1, 0);
  int iend = min(i0 + 1, this->nx - 1);
  int jstart = max(j0 - 1, 0);
  int jend = min(j0 + 1, this->ny - 1);
  int kstart = max(k0 - 1, 0);
  int kend = min(k0 + 1, this->nz - 1);

  // Loop over neighbouring bins
  for (int i = istart; i <= iend; i++) {
    for (int j = jstart; j <= jend; j++) {
      for (int k = kstart; k <= kend; k++) {

        // For each atom in the current bin, calculate the actual distance
        vector<int> binIndices = this->bins[i][j][k];
        for (auto &idx : binIndices) {
          double ix = this->positions(idx, 0);
          double iy = this->positions(idx, 1);
          double iz = this->positions(idx, 2);
          double deltax = x - ix;
          double deltay = y - iy;
          double deltaz = z - iz;
          double distanceSquared =
              deltax * deltax + deltay * deltay + deltaz * deltaz;
          if (distanceSquared <= this->cutoffSquared) {
            neighbours.push_back(idx);
            distancesSquared.push_back(distanceSquared);
            distances.push_back(sqrt(distanceSquared));
          }
        }
      }
    }
  }
  return CellListResult{neighbours, distances, distancesSquared};
}

CellListResultComplete
CellList::getAllNeighboursInfoForPosition(const double x, const double y,
                                          const double z) const {
  // The indices of the neighbouring atoms
  CellListResultComplete toreturn;
  constexpr size_t reserveGuess = 100;
  toreturn.indices.reserve(reserveGuess);
  toreturn.distancesSquared.reserve(reserveGuess);
  toreturn.distances.reserve(reserveGuess);
  toreturn.dx.reserve(reserveGuess);
  toreturn.dy.reserve(reserveGuess);
  toreturn.dz.reserve(reserveGuess);
  toreturn.dxSquared.reserve(reserveGuess);
  toreturn.dySquared.reserve(reserveGuess);
  toreturn.dzSquared.reserve(reserveGuess);
  // Find bin for the given position
  int i0 = (x - this->xmin) / this->dx;
  int j0 = (y - this->ymin) / this->dy;
  int k0 = (z - this->zmin) / this->dz;

  // Get the bin ranges to check for each dimension.
  int istart = max(i0 - 1, 0);
  int iend = min(i0 + 1, this->nx - 1);
  int jstart = max(j0 - 1, 0);
  int jend = min(j0 + 1, this->ny - 1);
  int kstart = max(k0 - 1, 0);
  int kend = min(k0 + 1, this->nz - 1);

  // Loop over neighbouring bins
  for (int i = istart; i <= iend; i++) {
    for (int j = jstart; j <= jend; j++) {
      for (int k = kstart; k <= kend; k++) {

        // For each atom in the current bin, calculate the actual distance
        for (auto &idx : this->bins[i][j][k]) {
          double deltax = this->positions(idx, 0) - x;
          double deltay = this->positions(idx, 1) - y;
          double deltaz = this->positions(idx, 2) - z;
          double distanceSquared =
              deltax * deltax + deltay * deltay + deltaz * deltaz;
          if (distanceSquared <= this->cutoffSquared) {
            toreturn.indices.emplace_back(idx);
            toreturn.distancesSquared.emplace_back(distanceSquared);
            toreturn.distances.emplace_back(sqrt(distanceSquared));
            toreturn.dx.emplace_back(deltax);
            toreturn.dy.emplace_back(deltay);
            toreturn.dz.emplace_back(deltaz);
            toreturn.dxSquared.emplace_back(deltax * deltax);
            toreturn.dySquared.emplace_back(deltay * deltay);
            toreturn.dzSquared.emplace_back(deltaz * deltaz);
          }
        }
      }
    }
  }
  return toreturn;
}

CellListResult CellList::getNeighboursForIndex(const int idx) const {
  double x = this->positions(idx, 0);
  double y = this->positions(idx, 1);
  double z = this->positions(idx, 2);
  CellListResult result = this->getNeighboursForPosition(x, y, z);

  // Remove self from neighbours
  for (int i = 0; i < result.indices.size(); ++i) {
    if (result.indices[i] == idx) {
      result.indices.erase(result.indices.begin() + i);
      result.distances.erase(result.distances.begin() + i);
      result.distancesSquared.erase(result.distancesSquared.begin() + i);
      break;
    }
  }
  return result;
}
