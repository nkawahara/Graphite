#include <stdlib.h>
#include <math.h>
//#include "network_node.h"
#include "network_model_auto_hop_counter.h"
#include "simulator.h"
#include "config.h"
#include "config.h"
#include "tile.h"
#include "constants.h"

void addEdge(int v, int u, int weight, struct NodeEdge *node){
        node[ u ].to.push_back( v );
        node[ u ].cost.push_back( weight );
        //有向グラフならここから下の処理が不要
        node[ v ].to.push_back( u );
        node[ v ].cost.push_back( weight );
}

void dijkstra(int n, int start, int end, struct NodeEdge *node){
        for(int i=0 ; i<n ; i++){
                node[i].done = false;
                node[i].minCost = -1;
        }

        node[start].minCost = 0;        //スタートノードまでのコストは0
        while(1){
                int donenode = -1;      //最新の確定したノード番号(-1はNULLのかわり)
                for(int i=0 ; i<n ; i++){

                        if( node[i].done==true ){//ノードiが確定しているときcontinue
                                continue;
                        }
                        if( node[i].minCost < 0 ){//ノードiまでの現時点での最小コストが不明のとき
                                continue;
                        }

                        //確定したノード番号が-1かノードiの現時点の最小コストが小さいとき
                        //確定ノード番号を更新する
                        if( donenode<0 || node[i].minCost < node[donenode].minCost){
                                donenode = i;
                        }
                }

                if(donenode==-1) break; //すべてのノードが確定したら終了

                node[donenode].done = true;     //ノードを確定させる

                for(unsigned int i=0 ; i<node[donenode].to.size() ; i++){
                        int to = node[donenode].to[i];
                        int cost = node[donenode].minCost + node[donenode].cost[i];

                        if( node[to].minCost < 0 || cost < node[to].minCost ){
                                node[to].minCost = cost;
                                node[to].from = donenode;

                       }
                }
        }
}


NetworkModelAutoHopCounter::NetworkModelAutoHopCounter(Network *net, SInt32 network_id)
   : NetworkModel(net, network_id)
   , _router_power_model(NULL)
   , _electrical_link_power_model(NULL)
{
_application_tiles = 16;
   assert(_application_tiles = Config::getSingleton()->getApplicationTiles());
   //_application_tiles = Config::getSingleton()->getApplicationTiles();
  
   try
   {
      _flit_width = Sim()->getCfg()->getInt("network/auto_hop_counter/flit_width");
   }
   catch (...)
   {
      LOG_PRINT_ERROR("Could not read emesh_hop_counter paramters from the cfg file");
   }

   // Broadcast Capability
   _has_broadcast_capability = false;

   createRouterAndLinkModels();
   
   // Initialize event counters
   initializeEventCounters();
}

NetworkModelAutoHopCounter::~NetworkModelAutoHopCounter()
{
   // Destroy the Router & Link Models
   destroyRouterAndLinkModels();
}

void
NetworkModelAutoHopCounter::createRouterAndLinkModels()
{
   if (isSystemTile(_tile_id))
      return;

   // Link parameters
   UInt64 link_delay = 0;
   string link_type;
   double link_length = 0.0;
   // Router parameters
   UInt64 router_delay = 0;   // Delay of the router (in clock cycles)
   UInt32 num_flits_per_output_buffer = 0;   // Here, contention is not modeled
    
   try
   {
      link_delay = (UInt64) Sim()->getCfg()->getInt("network/auto_hop_counter/link/delay");
      link_type = Sim()->getCfg()->getString("network/auto_hop_counter/link/type");
      link_length = Sim()->getCfg()->getFloat("general/tile_width");

      router_delay = (UInt64) Sim()->getCfg()->getInt("network/auto_hop_counter/router/delay");
      num_flits_per_output_buffer = Sim()->getCfg()->getInt("network/auto_hop_counter/router/num_flits_per_port_buffer");
   }
   catch (...)
   {
      LOG_PRINT_ERROR("Could not read emesh_hop_counter link and router parameters");
   }

   LOG_ASSERT_ERROR(link_delay == 1, "Network Link Delay(%llu) is not 1 cycle", link_delay);

   // Hop latency
   _hop_latency = router_delay + link_delay;
   
   // Router & Link have the same throughput (flit_width = phit_width = link_width)
   // Router & Link are clocked at the same frequency

   // define network.
   _route_flag = true;
   int cost = link_delay;
   addEdge(0,1,cost,_node);
    
   // Swiches Parametor
   // _sw = 45324421;

   //offset
   _offsetA = 1;
_sw = 2;
addEdge(1,17,cost,_node);
addEdge(2,17,cost,_node);
addEdge(3,17,cost,_node);
addEdge(4,17,cost,_node);
addEdge(5,17,cost,_node);
addEdge(6,17,cost,_node);
addEdge(7,17,cost,_node);
addEdge(8,17,cost,_node);
addEdge(9,18,cost,_node);
addEdge(10,18,cost,_node);
addEdge(11,18,cost,_node);
addEdge(12,18,cost,_node);
addEdge(13,18,cost,_node);
addEdge(14,18,cost,_node);
addEdge(15,18,cost,_node);
addEdge(16,18,cost,_node);
addEdge(17,18,cost,_node);
   // Instantiate router & link power models
   UInt32 num_router_ports = 5;
   if (Config::getSingleton()->getEnablePowerModeling())
   {
      _router_power_model = new RouterPowerModel(_frequency, _voltage, num_router_ports, num_router_ports,
                                                 num_flits_per_output_buffer, _flit_width);
      _electrical_link_power_model = new ElectricalLinkPowerModel(link_type, _frequency, _voltage, link_length, _flit_width);
   }

}

void
NetworkModelAutoHopCounter::destroyRouterAndLinkModels()
{
   if (isSystemTile(_tile_id))
      return;
   
   if (Config::getSingleton()->getEnablePowerModeling())
   {
      delete _router_power_model;
      delete _electrical_link_power_model;
   }
}

void
NetworkModelAutoHopCounter::initializeEventCounters()
{
   _buffer_writes = 0;
   _buffer_reads = 0;
   _switch_allocator_traversals = 0;
   _crossbar_traversals = 0;
   _link_traversals = 0;
}

void
NetworkModelAutoHopCounter::updateEventCounters(UInt32 num_flits, UInt32 num_hops)
{
   _buffer_writes += (num_flits * num_hops); 
   _buffer_reads += (num_flits * num_hops); 
   _switch_allocator_traversals += num_hops;
   _crossbar_traversals += (num_flits * num_hops);
   _link_traversals += (num_flits * num_hops);
}

void
NetworkModelAutoHopCounter::computePosition(tile_id_t tile, SInt32 &x, SInt32 &y)
{
   x = tile % _mesh_width;
   y = tile / _mesh_width;
}

SInt32
NetworkModelAutoHopCounter::computeDistance(SInt32 x1, SInt32 y1, SInt32 x2, SInt32 y2)
{
   return abs(x1 - x2) + abs(y1 - y2);
}

void
NetworkModelAutoHopCounter::routePacket(const NetPacket &pkt, queue<Hop> &next_hops)
{
   int start = (int)TILE_ID(pkt.sender);
   int end = (int)TILE_ID(pkt.receiver);
   UInt32 num_hops;
   //SInt32 path[];   
   vector<int> path;

   if ( _route_flag ) {
     // dijkstra( _application_tiles, start , end , _node );
    dijkstra( _application_tiles + _sw , start , end , _node );

    _route_flag=false;
   }
                                                                                                                           
   for(int i = end ; i != start ; i = _node[i].from ){
     path.push_back(i);
   }
   path.push_back(start);
   
   //最短経路の出力 
   cout << "最短経路は" << endl; 
   for(int i = path.size()-1 ; i >= 0 ; i--){
     cout << path[i] << " ";
   }
   cout << endl; 
   
   num_hops =  _node[end].minCost;
   //printf("%d -> %d =%d\n", start , end, num_hops);
   Latency latency = (isModelEnabled(pkt)) ? Latency(num_hops * _hop_latency,_frequency) : Latency(0,_frequency);

   updateDynamicEnergy(pkt, num_hops);

   Hop hop(pkt, TILE_ID(pkt.receiver), RECEIVE_TILE, latency, Time(0));
   next_hops.push(hop);
}

void
NetworkModelAutoHopCounter::outputSummary(std::ostream &out, const Time& target_completion_time)
{
   NetworkModel::outputSummary(out, target_completion_time);
   outputPowerSummary(out, target_completion_time);
   outputEventCountSummary(out);
}

// Power/Energy related functions
void
NetworkModelAutoHopCounter::updateDynamicEnergy(const NetPacket& packet, UInt32 num_hops)
{
   if (!isModelEnabled(packet))
      return;

   UInt32 num_flits = computeNumFlits(getModeledLength(packet));
    
   // Update event counters 
   updateEventCounters(num_flits, num_hops);

   // Update energy counters
   if (Config::getSingleton()->getEnablePowerModeling())
   {
      _router_power_model->updateDynamicEnergy(num_flits*num_hops, num_hops);
      _electrical_link_power_model->updateDynamicEnergy(num_flits * num_hops);
   }
}

void
NetworkModelAutoHopCounter::outputPowerSummary(ostream& out, const Time& target_completion_time)
{
   if (!Config::getSingleton()->getEnablePowerModeling())
      return;

   out << "    Power Model Statistics: " << endl;
   if (isApplicationTile(_tile_id))
   {
      // Convert time into seconds
      double target_completion_sec = target_completion_time.toSec();
      
      // Compute the final leakage/dynamic energy
      computeEnergy(target_completion_time);
      
      // We need to get the power of the router + all the outgoing links (a total of 4 outputs)
      double static_energy = _router_power_model->getStaticEnergy() +
                            (_electrical_link_power_model->getStaticEnergy() * _NUM_OUTPUT_DIRECTIONS);
      double dynamic_energy = _router_power_model->getDynamicEnergy() +
                              _electrical_link_power_model->getDynamicEnergy();
      out << "      Average Static Power (in W): " << static_energy / target_completion_sec << endl;
      out << "      Average Dynamic Power (in W): " << dynamic_energy / target_completion_sec << endl;
      out << "      Total Static Energy (in J): " << static_energy << endl;
      out << "      Total Dynamic Energy (in J): " << dynamic_energy << endl;
   }
   else if (isSystemTile(_tile_id))
   {
      out << "      Average Static Power (in W): " << endl;
      out << "      Average Dynamic Power (in W): " << endl;
      out << "      Static Energy (in J): " << endl;
      out << "      Dynamic Energy (in J): " << endl;
   }
   else
   {
      LOG_PRINT_ERROR("Unrecognized Tile ID(%i)", _tile_id);
   }
}

void
NetworkModelAutoHopCounter::outputEventCountSummary(ostream& out)
{
   out << "    Event Counters:" << endl;
   if (isApplicationTile(_tile_id))
   {
      out << "      Buffer Writes: " << _buffer_writes << endl;
      out << "      Buffer Reads: " << _buffer_reads << endl;
      out << "      Switch Allocator Traversals: " << _switch_allocator_traversals << endl;
      out << "      Crossbar Traversals: " << _crossbar_traversals << endl;
      out << "      Link Traversals: " << _link_traversals << endl;
   }
   else if (isSystemTile(_tile_id))
   {
      out << "      Buffer Writes: " << endl;
      out << "      Buffer Reads: " << endl;
      out << "      Switch Allocator Traversals: " << endl;
      out << "      Crossbar Traversals: " << endl;
      out << "      Link Traversals: " << endl;
   }
   else
   {
      LOG_PRINT_ERROR("Unrecognized Tile ID(%i)", _tile_id);
   }
}

void
NetworkModelAutoHopCounter::setDVFS(double frequency, double voltage, const Time& curr_time)
{
   if (!Config::getSingleton()->getEnablePowerModeling())
      return;

   LOG_PRINT("setDVFS[Frequency(%g), Voltage(%g), Time(%llu ns)] begin", frequency, voltage, curr_time.toNanosec());
   _router_power_model->setDVFS(frequency, voltage, curr_time);
   _electrical_link_power_model->setDVFS(frequency, voltage, curr_time);
   LOG_PRINT("setDVFS[Frequency(%g), Voltage(%g), Time(%llu ns)] end", frequency, voltage, curr_time.toNanosec());
}

void
NetworkModelAutoHopCounter::computeEnergy(const Time& curr_time)
{
   assert (Config::getSingleton()->getEnablePowerModeling());
   _router_power_model->computeEnergy(curr_time);
   _electrical_link_power_model->computeEnergy(curr_time);
}

double
NetworkModelAutoHopCounter::getDynamicEnergy()
{
   assert (Config::getSingleton()->getEnablePowerModeling());
   double dynamic_energy = _router_power_model->getDynamicEnergy() +
                           _electrical_link_power_model->getDynamicEnergy();
   return dynamic_energy;
}

double
NetworkModelAutoHopCounter::getStaticEnergy()
{
   assert (Config::getSingleton()->getEnablePowerModeling());
   double static_energy = _router_power_model->getStaticEnergy() +
                          (_electrical_link_power_model->getStaticEnergy() * _NUM_OUTPUT_DIRECTIONS);
   return static_energy;
}
