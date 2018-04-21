/******************************************************************************
 * Copyright (c) 2016-2017, Bradley J Chambers (brad.chambers@gmail.com)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of Hobu, Inc. or Flaxen Geo Consulting nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include "DistanceToLineFilter.hpp"

#include <pdal/EigenUtils.hpp>
#include <pdal/KDIndex.hpp>
#include <pdal/util/ProgramArgs.hpp>

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace pdal
{

static PluginInfo const s_info =
    PluginInfo("filters.distance_to_polyline_filter", "Distance To Polyline Filter",
               "http://pdal.io/stages/filters.distance_to_polyline.html");

CREATE_STATIC_PLUGIN(1, 0, DistanceToLineFilter, Filter, s_info)


std::string DistanceToLineFilter::getName() const
{
    return s_info.name;
}

void DistanceToLineFilter::addArgs(ProgramArgs& args)
{
    Arg& candidate = args.add("candidate", "candidate file name",
        m_candidateFile);
    args.add("layer", "Datasource layer to use", m_layer);
}

void DistanceToLineFilter::addDimensions(PointLayoutPtr layout)
{
    using namespace Dimension;

    layout->registerDims(
        {Id::LineDirX, Id::LineDirY, Id::LineDirZ, Id::LineDist});
}

void DistanceToLineFilter::prepared(PointTableRef table)
{
}

std::vector<double> DistanceToLineFilter::calcDisplacement(PointRef p1, PointRef p2)
{
    std::vector<double> ret;
    ret.push_back(p1.getFieldAs<double>(Dimension::Id::X) - p2.getFieldAs<double>(Dimension::Id::X));
    ret.push_back(p1.getFieldAs<double>(Dimension::Id::Y) - p2.getFieldAs<double>(Dimension::Id::Y));
    ret.push_back(p1.getFieldAs<double>(Dimension::Id::Z) - p2.getFieldAs<double>(Dimension::Id::Z));
    return ret;
}

PointViewPtr DistanceToLineFilter::loadSet(const std::string& filename,
    PointTable& table)
{
    PipelineManager mgr;

    Stage& candidateReader = m_candReader ? *m_candReader : mgr.makeReader(filename, "");
    candidateReader.prepare(table);
    PointViewSet viewSet = candidateReader.execute(table);
    assert(viewSet.size() == 1);
    return *viewSet.begin();
}

void DistanceToLineFilter::filter(PointView& view)
{
    PointTable candTable;
    PointViewPtr candView = loadSet(m_candidateFile, candTable);
    KD3Index kdiCand(*candView);
    kdiCand.build();

    PointRef point_src(view, 0);
    for (PointId i = 0; i < view.size(); ++i)
    {
        // Find the closest two points on the line.
        point_src.setPointId(i);
        std::vector<PointId> nn_id = kdiCand.neighbors(point_src, 2);

        // Calculate the displacement vector to the line segment
        Eigen::Vector3f A(view.getFieldAs<double>(Dimension::Id::X, i),
          view.getFieldAs<double>(Dimension::Id::Y, i),
          view.getFieldAs<double>(Dimension::Id::Z, i));
        Eigen::Vector3f B(candView->getFieldAs<double>(Dimension::Id::X, nn_id[0]), 
          candView->getFieldAs<double>(Dimension::Id::Y, nn_id[0]),
          candView->getFieldAs<double>(Dimension::Id::Z, nn_id[0]));
        Eigen::Vector3f C(candView->getFieldAs<double>(Dimension::Id::X, nn_id[1]),
          candView->getFieldAs<double>(Dimension::Id::Y, nn_id[1]),
          candView->getFieldAs<double>(Dimension::Id::Z, nn_id[1]));
        Eigen::Vector3f u = eigen::computeDisplacement(A, B, C);
        //std::cout << "DistanceToLineFilter::filter().\n";
        //std::cout << "  Source point is " << A.transpose() << "\n";
        //std::cout << "  Two nearest points are " << B.transpose() << " and " << C.transpose() << "\n";
        //std::cout << "  Displacement Vector is " << u.transpose() << "\n";

        view.setField(Dimension::Id::LineDirX, i, u[0]);
        view.setField(Dimension::Id::LineDirY, i, u[1]);
        view.setField(Dimension::Id::LineDirZ, i, u[2]);

        float dist = u.norm();
        view.setField(Dimension::Id::LineDist, i, dist);
    }
}

} // namespace pdal
