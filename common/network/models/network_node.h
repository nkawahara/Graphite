#include <vector>

struct NodeEdge{
        vector<int> to;         //どのノードとつながっているか
        vector<int> cost;       //エッジのコスト

        //ここから下はダイクストラ法のために必要な情報
        bool done;              //確定ノードかどうか
        int minCost;    //スタートノードからのコスト
        int from;               //どのノードから来たか
};
