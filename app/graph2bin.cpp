/******************************************************************************
 * graph2bin.cpp
 *
 * Converts a graph in the METIS text format into the binary form that
 * graph_io::readGraphBinary reads.
 *
 * Reading the text format is dominated by turning characters into integers,
 * not by the disk: on a 3.6 GB instance the bytes arrive in half a second and
 * the parsing takes twenty-three. The binary form is the same graph as a CSR
 * dump, so loading it is a few reads and no parsing.
 *
 *   graph2bin <input.metis> <output.bin>
 *
 * The output is byte for byte the graph the reader built, so a run on the
 * converted file sees exactly the graph a run on the text file would.
 *****************************************************************************/

#include <iostream>
#include <memory>
#include <string>

#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    if (argn != 3) {
        std::cerr << "usage: " << argv[0]
                  << " <input.metis> <output.bin>" << std::endl;
        return 2;
    }
    const std::string in = argv[1];
    const std::string out = argv[2];

    timer t;
    auto G = graph_io::readGraphWeighted<mutable_graph>(in);
    const double t_read = t.elapsed();
    if (!G) {
        std::cerr << "could not read " << in << std::endl;
        return 1;
    }

    t.restart();
    int rc = graph_io::writeGraphBinary<mutable_graph>(G, out);
    const double t_write = t.elapsed();
    if (rc != 0) {
        std::cerr << "could not write " << out << std::endl;
        return 1;
    }

    std::cout << "converted " << in
              << " n=" << G->number_of_nodes()
              << " m=" << G->number_of_edges() / 2
              << " read=" << t_read << "s"
              << " write=" << t_write << "s" << std::endl;
    return 0;
}
