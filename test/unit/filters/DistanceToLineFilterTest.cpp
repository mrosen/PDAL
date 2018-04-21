/******************************************************************************
* Copyright (c) 2015, Hobu Inc. (hobu@hobu.co)
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

#include <pdal/pdal_test_main.hpp>

#include <io/FauxReader.hpp>
#include <filters/DistanceToLineFilter.hpp>

using namespace pdal;

TEST(DistanceToLineFilterTest, easy_debugging)
{
    // This isn't a real test.  It's just here to allow easy debugging.

    point_count_t count = 3;

    Options readerOps;
    readerOps.add("bounds", BOX3D(1, 1, 1, count, count, count));
    readerOps.add("mode", "ramp");
    readerOps.add("count", count);

    //FauxReader r;
    std::shared_ptr<FauxReader> r(new FauxReader);
    r->setOptions(readerOps);

    PointTable t;
    r->prepare(t);
    PointViewSet s = r->execute(t);

    EXPECT_EQ(s.size(), 1u);
    PointViewPtr v = *s.begin();
    EXPECT_EQ(v->size(), (size_t)count);

    //for (PointId i = 0; i < count; i++)
    //    std::cerr << "X[" << i << "] = " <<
    //        v->getFieldAs<double>(Dimension::Id::X, i) << "!\n";
}

TEST(DistanceToLineFilterTest, trivial)
{
    // Smallest possible test:  CL 0-1 and a PC of a single point

    Options rOps;
    rOps.replace("bounds", BOX3D(0, 0, 0, 1, 0, 0));
    rOps.replace("mode", "ramp");
    rOps.replace("count", 2);
    //FauxReader rCenterLine;
    std::shared_ptr<FauxReader> rCenterLine(new FauxReader);
    rCenterLine->setOptions(rOps);

    rOps.replace("bounds", BOX3D(.25, 1, 0, .25, 1, 0));
    rOps.replace("mode", "ramp");
    rOps.replace("count", 1);

    FauxReader rPointCloud;
    rPointCloud.setOptions(rOps);

    DistanceToLineFilter f;
    f.setInput(rPointCloud);
    f.setCandidateReader(rCenterLine);

    PointTable t;
    f.prepare(t);
    PointViewSet s = f.execute(t);

    point_count_t count = 1;
    EXPECT_EQ(s.size(), 1u);
    PointViewPtr v = *s.begin();
    EXPECT_EQ(v->size(), (size_t)count);

    //for (PointId i = 0; i < count; i++)
    //    std::cerr << "X[" << i << "] = " <<
    //        v->getFieldAs<double>(Dimension::Id::X, i) << "!\n";

    EXPECT_EQ(v->getFieldAs<double>(Dimension::Id::LineDirX, 0), 0.0);
    EXPECT_EQ(v->getFieldAs<double>(Dimension::Id::LineDirY, 0), -1.0);
    EXPECT_EQ(v->getFieldAs<double>(Dimension::Id::LineDirZ, 0), 0.0);
}
void ConfigureFilterForTest(const BOX3D &boxPC, const BOX3D &boxCL, size_t nPC, size_t nCL, DistanceToLineFilter &f, Options &rOps, FauxReader *rPointCloud)
{

    // PointCloud
    rOps.replace("bounds", boxPC);
    rOps.replace("mode", "ramp");
    rOps.replace("count", nPC);
    rPointCloud->setOptions(rOps);

    // Vector Centerline
    rOps.replace("bounds", boxCL);
    rOps.replace("mode", "ramp");
    rOps.replace("count", nCL);
    std::shared_ptr<FauxReader> rCenterLine(new FauxReader);
    rCenterLine->setOptions(rOps);

    f.setInput(*rPointCloud);
    f.setCandidateReader(rCenterLine);
}

TEST(DistanceToLineFilterTest, simple)
{
    // Simple example lifted from 
    // https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d

    DistanceToLineFilter f;
    Options rOps;
    FauxReader rPointCloud;
    ConfigureFilterForTest(BOX3D(4.0, 2.0, 1.0, 4.0, 2.0, 1.0), BOX3D(1, 0, 1, 1, 2, 0), 1, 2, f, rOps, &rPointCloud); 

    PointTable t;
    f.prepare(t);
    PointViewSet s = f.execute(t);

    PointViewPtr v = *s.begin();
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirX, 0), -3.0);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirY, 0), -0.4);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirZ, 0), -0.8);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDist, 0), 3.1304953);
}
TEST(DistanceToLineFilterTest, proj_not_avail)
{
    // Edge case: can't project A onto BC b/c BC is too short.

    DistanceToLineFilter f;
    Options rOps;
    FauxReader rPointCloud;
    // A is at X = 1.5; BC runs from X = 0 to 1; Projection is X = 1.5 but
    // closest point is 1.0
    ConfigureFilterForTest(BOX3D(1.5, 0, 0, 1.5, 0, 0), BOX3D(0, 1, 0, 1, 1, 0), 1, 2, f, rOps, &rPointCloud); 

    PointTable t;
    f.prepare(t);
    PointViewSet s = f.execute(t);

    PointViewPtr v = *s.begin();
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirX, 0), -0.5);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirY, 0), 1);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDirZ, 0), 0);
    EXPECT_FLOAT_EQ(v->getFieldAs<float>(Dimension::Id::LineDist, 0), 1.118034);
}
