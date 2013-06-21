/* 

Example of orthogonal embeddings via OGDF.  The trick is you need to set the
node size or else it will through a cryptic "floating exception".  See 

http://www.ogdf.net/forum/showthread.php?tid=427

for details.  

*/

/*
  Various other simple layouts that work
  
  #include <ogdf/energybased/FMMMLayout.h>
  FMMMLayout fmmm;
  fmmm.call(GA);
  
  #include <ogdf/planarlayout/SchnyderLayout.h>
  SchnyderLayout sl;
  sl.call(GA);
*/


#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/misclayout/CircularLayout.h>

using namespace ogdf;
 




int main()
{
    Graph G;
    // planarTriconnectedGraph(G, 200, 300)
    GraphAttributes GA(G);
    GA.readGML(G, "graphs/link.gml");
    node v;
    forall_nodes(v, G) {
	if (v->degree() == 2){
	    GA.width(v) = 1.0;
	    GA.height(v) = 1.0;
	}
	else{
	    GA.width(v) = 1.0;
	    GA.height(v) = 1.0;
	}
    }

    PlanarizationLayout pl;
    OrthoLayout *ol = new OrthoLayout;
    pl.setPlanarLayouter(ol);
    pl.call(GA);

    // GA.writeGML("graphs/er-diagram-layout2.gml");
    GA.writeSVG("graphs/link.svg");
    return 0;
}
