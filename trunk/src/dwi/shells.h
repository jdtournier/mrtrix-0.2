/*
    Copyright 2013 Brain Research Institute, Melbourne, Australia

    Written by B Jeurissen, 12/08/13.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __dwi_shells_h__
#define __dwi_shells_h__

#include "math/matrix.h"
#include <list>
#include <limits>

namespace MR
{
  namespace DWI
  {

    class Shell
    {

    public:
      Shell()
      {
        count_ = 0;
        avg_bval_ = 0;
        std_bval_ = 0;
        min_bval_ = 0;
        max_bval_ = 0;
      }

      Shell(std::vector<double> bvals, std::vector<int> clusters, int cluster)
      {
        count_ = 0;
        avg_bval_ = 0;
        std_bval_ = 0;
        min_bval_ = std::numeric_limits<double>::max();
        max_bval_ = std::numeric_limits<double>::min();
        for (size_t j = 0; j < clusters.size(); j++) {
          if (clusters[j] == cluster) {
            idx_.push_back(j);
            count_++;
            avg_bval_ += bvals[j];
            if (bvals[j]<min_bval_) {
              min_bval_ = bvals[j];
            }
            if (bvals[j]>max_bval_) {
              max_bval_ = bvals[j];
            }
          }
        }
        avg_bval_ = avg_bval_/count_;
        for (size_t j = 0; j < clusters.size(); j++) {
          if (clusters[j] == cluster) {
            std_bval_ += pow(bvals[j]-avg_bval_,2);
          }
        }
        std_bval_ = sqrt(std_bval_/count_);
      }

      std::vector<int> idx()
      {
        return idx_;
      }

      int count() const
      {
        return count_;
      }

      double avg_bval() const
      {
        return avg_bval_;
      }

      double std_bval() const
      {
        return std_bval_;
      }

      double min_bval() const
      {
        return min_bval_;
      }

      double max_bval() const
      {
        return max_bval_;
      }

      Shell& operator= (const Shell& rhs)
      {
        idx_ = rhs.idx_;
        avg_bval_ = rhs.avg_bval_;
        std_bval_ = rhs.std_bval_;
        count_ = rhs.count_;
        min_bval_ = rhs.min_bval_;
        max_bval_ = rhs.max_bval_;
        return *this;
      }

      friend std::ostream& operator <<(std::ostream& stream, const Shell& S)
      {
        stream << S.count() << std::endl;
        stream << S.min_bval() << std::endl;
        stream << S.max_bval() << std::endl;
        stream << S.avg_bval() << std::endl;
        stream << S.std_bval() << std::endl;
        return stream;
      }

    private:
      std::vector<int> idx_;
      double avg_bval_;
      double std_bval_;
      int count_;
      double min_bval_;
      double max_bval_;
    };


    class Shells
    {
    public:
      Shells(const Math::Matrix& grad, size_t minDirections = 6, double bvalue_threshold = NAN)
      {
        for (size_t i = 0; i < grad.rows(); i++) {
          bvals.push_back(grad (i,3));
        }
        minBval = *std::min_element(bvals.begin(),bvals.end());
        maxBval = *std::max_element(bvals.begin(),bvals.end());
        if (!finite (bvalue_threshold)) {
          bvalue_threshold = 100;
          //bvalue_threshold = (maxBval-minBval)/(2.0*double(bvals.size()));
          //if (!bvalue_threshold)
          //  bvalue_threshold = 1;
        }
        clusterBvalues(minDirections, bvalue_threshold);
        sortByBval();
      }

      int count()
      {
        return shells.size();
      }

      Shell& operator[] (const int i)
      {
        return shells[i];
      }

      Shell& first()
      {
        return shells.front();
      }

      Shell& last()
      {
        return shells.back();
      }

      void sortByCount()
      {
        std::sort (shells.begin(), shells.end(), countComp);
      }

      void sortByBval()
      {
        std::sort (shells.begin(), shells.end(), bvalComp);
      }

      friend std::ostream& operator <<(std::ostream& stream, const Shells& S)
      {
        for(typename std::vector<Shell>::const_iterator it = S.shells.begin(); it != S.shells.end(); ++it) {
          stream << *it << std::endl;
        }
        return stream;
      }

    private:
      std::vector<Shell> shells;
      std::vector<double> bvals;
      double minBval;
      double maxBval;

      static bool countComp (Shell a, Shell b)
      {
        return (a.count()<b.count());
      }

      static bool bvalComp (Shell a, Shell b)
      {
        return (a.avg_bval()<b.avg_bval());
      }

      void regionQuery (double p, std::vector<double> x, double eps, std::vector<int>& idx)
      {
        for (size_t i = 0; i < x.size(); i++) {
          if (std::abs(p-x[i]) < eps) {
            idx.push_back(i);
          }
        }
      }

      void clusterBvalues(size_t minDirections, double eps)
      {
        std::vector<bool> visited;
        visited.resize(bvals.size(),false);

        std::vector<int> cluster;
        cluster.resize(bvals.size(),-1);
        int clusterIdx = -1;

        for (size_t ii = 0; ii < bvals.size(); ii++) {
          if (!visited[ii]) {
            visited[ii] = true;
            std::vector<int> neighborIdx;
            regionQuery(bvals[ii],bvals,eps,neighborIdx);
            if (bvals[ii] > eps && neighborIdx.size() < minDirections) {
              cluster[ii] = -1;
            } else {
              cluster[ii]=++clusterIdx;
              for (size_t i = 0; i < neighborIdx.size(); i++) {
                if (!visited[neighborIdx[i]]) {
                  visited[neighborIdx[i]] = true;
                  std::vector<int> neighborIdx2;
                  regionQuery(bvals[neighborIdx[i]],bvals,eps,neighborIdx2);
                  if (neighborIdx2.size() >= minDirections) {
                    for (size_t j = 0; j < neighborIdx2.size(); j++) {
                      neighborIdx.push_back(neighborIdx2[j]);
                    }
                  }
                }
                if (cluster[neighborIdx[i]] < 0) {
                  cluster[neighborIdx[i]] = clusterIdx;
                }
              }
            }
          }
        }
        double minIdx = *std::min_element(cluster.begin(),cluster.end());
        double maxIdx = *std::max_element(cluster.begin(),cluster.end());
        for (size_t i = minIdx; i <= maxIdx; i++) {
          Shell s(bvals,cluster,i);
          shells.push_back(s);
        }
      }
    };
  }
}

#endif

