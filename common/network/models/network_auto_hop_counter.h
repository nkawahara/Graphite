#pragma once

#include "network.h"
#include "network_model.h"
#include "router_power_model.h"
#include "electrical_link_power_model.h"
#include "network_node.h"
#include "lock.h"
#include <vector>
#include <queue>

class NetworkModelAutoHopCounter : public NetworkModel
{
public:
   NetworkModelAutoHopCounter(Network *net, SInt32 network_id);
   ~NetworkModelAutoHopCounter();

   void routePacket(const NetPacket &pkt, queue<Hop> &next_hops);
   void outputSummary(std::ostream &out, const Time& target_completion_time);

   // Energy computation
   void computeEnergy(const Time& curr_time);
   double getDynamicEnergy();
   double getStaticEnergy();

private:
 //  SInt32 _edge_num = 20;
   // Topolgy route data
   struct NodeEdge _node[20];
   bool _route_flag;
   SInt32 _application_tiles;

   // Topolgy parameters
   SInt32 _mesh_width;
   SInt32 _mesh_height;
   static const UInt32 _NUM_OUTPUT_DIRECTIONS = 5;

   // 
   // network
   // 
   SInt32 _sw;
   SInt32 _offsetA;

   // Electrical router and link power models
   RouterPowerModel* _router_power_model;
   ElectricalLinkPowerModel* _electrical_link_power_model;
   // Latency parameters
   UInt64 _hop_latency;

   // Event counters
   UInt64 _buffer_writes;
   UInt64 _buffer_reads;
   UInt64 _switch_allocator_traversals;
   UInt64 _crossbar_traversals;
   UInt64 _link_traversals;

   // DVFS
   void setDVFS(double frequency, double voltage, const Time& curr_time);
   
   // Create/destroy router/link models
   void createRouterAndLinkModels();
   void initializeEventCounters();
   void destroyRouterAndLinkModels();
   
   void computePosition(tile_id_t tile, SInt32 &x, SInt32 &y);
   SInt32 computeDistance(SInt32 x1, SInt32 y1, SInt32 x2, SInt32 y2);
   void updateDynamicEnergy(const NetPacket& packet, UInt32 num_hops);
   void updateEventCounters(UInt32 num_flits, UInt32 num_hops);
   
   // Summary
   void outputPowerSummary(ostream& out, const Time& target_completion_time);
   void outputEventCountSummary(ostream& out);

   // add_on dijikstra
   //UInt64 dijkstra(int s, int flag);
};
