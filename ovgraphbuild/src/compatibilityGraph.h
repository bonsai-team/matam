#ifndef COMPATIBILTYGRAPH_H
#define COMPATIBILTYGRAPH_H

#include <seqan/graph_types.h>
//#include <seqan/graph_algorithms.h>

using TGraph = seqan::Graph<seqan::Directed<> >;
using TVertexDescriptor = seqan::VertexDescriptor<TGraph>::Type;
using TEdgeDescriptor = seqan::EdgeDescriptor<TGraph>::Type;
using TSize = seqan::Size<TGraph>::Type;

using TProperties = seqan::String<seqan::CharString>;

#endif // COMPATIBILTYGRAPH_H
