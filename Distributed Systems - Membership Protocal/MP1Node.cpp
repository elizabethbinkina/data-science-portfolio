/**********************************
 * FILE NAME: MP1Node.cpp
 *
 * DESCRIPTION: Membership protocol run by this Node.
 * 				Definition of MP1Node class functions. (Revised 2020)
 *
 *  Starter code template
 **********************************/

#include "MP1Node.h"
#include "Member.h"

//ADDED 

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <utility>

/*
 * Note: You can change/add any functions in MP1Node.{h,cpp}
 */

/**
 * Overloaded Constructor of the MP1Node class
 * You can add new members to the class if you think it
 * is necessary for your logic to work
 */
MP1Node::MP1Node( Params *params, EmulNet *emul, Log *log, Address *address) {
	for( int i = 0; i < 6; i++ ) {
		NULLADDR[i] = 0;
	}
	this->memberNode = new Member;
    this->shouldDeleteMember = true;
	memberNode->inited = false;
	this->emulNet = emul;
	this->log = log;
	this->par = params;
	this->memberNode->addr = *address;
}

/**
 * Overloaded Constructor of the MP1Node class
 * You can add new members to the class if you think it
 * is necessary for your logic to work
 */
MP1Node::MP1Node(Member* member, Params *params, EmulNet *emul, Log *log, Address *address) {
    for( int i = 0; i < 6; i++ ) {
        NULLADDR[i] = 0;
    }
    this->memberNode = member;
    this->shouldDeleteMember = false;
    this->emulNet = emul;
    this->log = log;
    this->par = params;
    this->memberNode->addr = *address;
}

/**
 * Destructor of the MP1Node class
 */
MP1Node::~MP1Node() {
    if (shouldDeleteMember) {
        delete this->memberNode;
    }
}

/**
* FUNCTION NAME: recvLoop
*
* DESCRIPTION: This function receives message from the network and pushes into the queue
*                 This function is called by a node to receive messages currently waiting for it
*/
int MP1Node::recvLoop() {
    if ( memberNode->bFailed ) {
    	return false;
    }
    else {
    	return emulNet->ENrecv(&(memberNode->addr), enqueueWrapper, NULL, 1, &(memberNode->mp1q));
    }
}

/**
 * FUNCTION NAME: enqueueWrapper
 *
 * DESCRIPTION: Enqueue the message from Emulnet into the queue
 */
int MP1Node::enqueueWrapper(void *env, char *buff, int size) {
	Queue q;
	return q.enqueue((queue<q_elt> *)env, (void *)buff, size);
}

/**
* FUNCTION NAME: nodeStart
*
* DESCRIPTION: This function bootstraps the node
*                 All initializations routines for a member.
*                 Called by the application layer.
*/
void MP1Node::nodeStart(char *servaddrstr, short servport) {
    Address joinaddr;
    joinaddr = getJoinAddress();

    // Self booting routines
    if( initThisNode(&joinaddr) == -1 ) {
#ifdef DEBUGLOG
        log->LOG(&memberNode->addr, "init_thisnode failed. Exit.");
#endif
        exit(1);
    }

    if( !introduceSelfToGroup(&joinaddr) ) {
        finishUpThisNode();
#ifdef DEBUGLOG
        log->LOG(&memberNode->addr, "Unable to join self to group. Exiting.");
#endif
        exit(1);
    }

    return;
}

/**
 * FUNCTION NAME: initThisNode
 *
 * DESCRIPTION: Find out who I am and start up
 */
int MP1Node::initThisNode(Address *joinaddr) {
    /*
    * This function is partially implemented and may require changes
    */

    //COMMENTED THIS OUT FOR NOW 

	//int id = *(int*)(&memberNode->addr.addr);
	//int port = *(short*)(&memberNode->addr.addr[4]);

	memberNode->bFailed = false;
	memberNode->inited = true;
	memberNode->inGroup = false;
	// node is up!
	memberNode->nnb = 0;
	memberNode->heartbeat = 0;

	memberNode->pingCounter = 0;
	memberNode->timeOutCounter = TFAIL;

    memberNode->pingCounter = 0;
    memberNode->timeOutCounter = TFAIL;

    initMemberListTable(memberNode);

    if (joinaddr != NULL) {
        //send JOINREQ to introducter
        size_t msgsize = sizeof(MessageHdr) + sizeof(Address) + sizeof(long); //addr + heartbeat THIS WAS the problem, needed to add sizeof(Adresss)
        MessageHdr *msg = (MessageHdr *) malloc(msgsize);
        msg->msgType = JOINREQ;

        char *payload = (char *)(msg + 1);
        memcpy(payload, &memberNode->addr, sizeof(Address));
        payload += sizeof(Address);
        memcpy(payload, &memberNode->heartbeat, sizeof(long));

        //memcpy(&msg->addr, &memberNode->addr, sizeof(memberNode->addr));
        //memcpy((char*)(msg+1), &memberNode->heartbeat, sizeof(long));

        emulNet->ENsend(&memberNode->addr, joinaddr, (char*)msg, msgsize);
        free(msg);
        
    } else {
        //means its the first node , add to self
        memberNode->inGroup = true;
        addSelfToMemberList();
    }
    return 0;
}


//Helper function for initThisNode

void MP1Node::addSelfToMemberList() {
    int id = 0;
    short port = 0;

    memcpy(&id, &memberNode->addr.addr[0], sizeof(int));
    memcpy(&port, &memberNode->addr.addr[4], sizeof(short));

    MemberListEntry selfEntry(id,port,memberNode->heartbeat, par->getcurrtime());
    memberNode->memberList.push_back(selfEntry);
    timeout_map[{id,port}] = 0;

}



/**
 * FUNCTION NAME: introduceSelfToGroup
 *
 * DESCRIPTION: Join the distributed system
 */
int MP1Node::introduceSelfToGroup(Address *joinaddr) {
	MessageHdr *msg;
#ifdef DEBUGLOG
    static char s[1024];
#endif

    if ( 0 == strcmp((char *)&(memberNode->addr.addr), (char *)&(joinaddr->addr))) {
        // I am the group booter (first process to join the group). Boot up the group
#ifdef DEBUGLOG
        log->LOG(&memberNode->addr, "Starting up group...");
#endif
        memberNode->inGroup = true;
    }
    else {
        size_t msgsize = sizeof(MessageHdr) + sizeof(joinaddr->addr) + sizeof(long);// + 1 ; //removed + 1
        msg = (MessageHdr *) malloc(msgsize * sizeof(char));

        // create JOINREQ message: format of data is {struct Address myaddr}
        msg->msgType = JOINREQ;
        memcpy((char *)(msg+1), &memberNode->addr.addr, sizeof(memberNode->addr.addr));
        memcpy((char *)(msg+1) + sizeof(memberNode->addr.addr), &memberNode->heartbeat, sizeof(long));

#ifdef DEBUGLOG
        sprintf(s, "Trying to join...");
        log->LOG(&memberNode->addr, s);
#endif

        // send JOINREQ message to introducer member
        emulNet->ENsend(&memberNode->addr, joinaddr, (char *)msg, msgsize);

        free(msg);
    }

    return 1;

}

/**
* FUNCTION NAME: finishUpThisNode
*
* DESCRIPTION: Wind up this node and clean up state
*/
int MP1Node::finishUpThisNode(){
    /*
     * Your code goes here
     */
     memberNode->inGroup = false;
     memberNode->memberList.clear();
     return 0 ;
}

/**
* FUNCTION NAME: nodeLoop
*
* DESCRIPTION: Executed periodically at each member
*                 Check your messages in queue and perform membership protocol duties
*/
void MP1Node::nodeLoop() {
    
    if (memberNode->bFailed) {
    	return;
    }

    // Check my messages
    checkMessages();

    // Wait until you're in the group...
    if( !memberNode->inGroup ) {
    	return;
    }

    // ...then jump in and share your responsibilites!
    nodeLoopOps();

    return;

}

/**
 * FUNCTION NAME: checkMessages
 *
 * DESCRIPTION: Check messages in the queue and call the respective message handler
 */
void MP1Node::checkMessages() {
    void *ptr;
    int size;

    // Pop waiting messages from memberNode's mp1q
    while ( !memberNode->mp1q.empty() ) {
    	ptr = memberNode->mp1q.front().elt;
    	size = memberNode->mp1q.front().size;
    	memberNode->mp1q.pop();
    	recvCallBack((void *)memberNode, (char *)ptr, size);
    }
    return;
}

/**
 * FUNCTION NAME: recvCallBack
 *
 * DESCRIPTION: Message handler for different message types
 */
bool MP1Node::recvCallBack(void *env, char *data, int size ) {   // check this
    /*
    * Your code goes here
    */

    MessageHdr *msg = (MessageHdr *) data;
    char *ptr = (char *)(msg + 1);


    switch (msg->msgType) {

        case JOINREQ: {
            Address senderAddr;
            memcpy(&senderAddr.addr, ptr, sizeof(Address));
            ptr += sizeof(Address);

            //extract heartbeat ---- THIS IS WHERE IM GOING WRONG
            long heartbeat;
            memcpy(&heartbeat, ptr, sizeof(long));  // THIS WAS THE PROBLEM (int64_t)

            //adding sender to membership lift if not present
            addOrUpdateMember(senderAddr, heartbeat, par->getcurrtime()); 
            
            //log additon
            //log->logNodeAdd(&memberNode->addr, &senderAddr); ---> removed for noe
            //send full membership list
            sendJoinReply(senderAddr);

            #ifdef DEBUGLOG
            log->LOG(&memberNode->addr, "DBG after JOINREQ:");
            debugPrintMemberList("after JOINREQ");
            #endif



            break;
        }

        case JOINREP: {

        
            //merge recived membership list
            mergeMembershipList(ptr, size - sizeof(MessageHdr));
            memberNode->inGroup = true;


            #ifdef DEBUGLOG
            log->LOG(&memberNode->addr, "DBG after JOINREP:");
            debugPrintMemberList("after JOINREP");
            #endif


            break;
        }

        case GOSSIP: {

            //merge received membership llist
            mergeMembershipList(ptr, size - sizeof(MessageHdr));

            #ifdef DEBUGLOG
            log->LOG(&memberNode->addr, "DBG after GOSSIP:");
            debugPrintMemberList("after GOSSIP");
            #endif


            break;
        }
        default:
            break;
    }

    return true;

}


// HELPER FUNCTIONS

void MP1Node::addOrUpdateMember(Address addr, long heartbeat, long timestamp) {
    int id;
    short port;
    memcpy(&id,&addr.addr[0], sizeof(int));
    memcpy(&port,&addr.addr[4], sizeof(short));

    bool found = false;

    for (auto &m : memberNode->memberList) {
        if(m.id == id && m.port == port) {
            if(heartbeat > m.heartbeat) {
                m.heartbeat = heartbeat;
                m.timestamp = timestamp; 
            }

            
            //reset timeout counter in the map
            timeout_map[{id, port}] = 0;
            found = true;
            break;
        }                      
        
    }

    if (!found) {
        MemberListEntry newEntry(id,port,heartbeat, timestamp);
        memberNode->memberList.push_back(newEntry);
        //initilize timeout counter in map 
        timeout_map[{id,port}] = 0;
        log->logNodeAdd(&memberNode->addr, &addr);
    }
}


void MP1Node::gossipMembershipList() {
    //choosing random neighbor
    int fanout = 3;
    //long currTime = par->getcurrtime();
    //vector<MemberListEntry> aliveList; // = memberNode->memberList;
    vector<MemberListEntry> allMembers = memberNode->memberList;
    random_shuffle(allMembers.begin(), allMembers.end());

    /*
    for (auto &m : memberNode->memberList) {
        if (currTime - m.timestamp <= TREMOVE) {
            aliveList.push_back(m);
        }
    }

    random_shuffle(aliveList.begin(), aliveList.end());
    */


    for (int i = 0; i < fanout && i < (int)allMembers.size(); i++) {
        Address target = getNodeAddress(allMembers[i].id, allMembers[i].port);

        //skip sending gossip to itself
        if(memcmp(&target.addr, &memberNode->addr.addr, 6) != 0) // THIS WAS THE PROBLEM (memcmp)
            sendGossip(target);
    }
}

// coverting id and port to adress
Address MP1Node::getNodeAddress(int id, short port) {
    Address addr;
    memcpy(&addr.addr[0], &id, sizeof(int));
    memcpy(&addr.addr[4], &port, sizeof(short));
    return addr;
}

//sending JOINREP w/ curr membership list to new node
void MP1Node::sendJoinReply(Address target) {
    //calc message size = header + membership list


    
    size_t listSize = memberNode->memberList.size() * sizeof(MemberListEntry);
    size_t msgSize = sizeof(MessageHdr) + listSize;

    MessageHdr *msg = (MessageHdr *)malloc(msgSize);
    msg->msgType = JOINREP;
    memcpy((char *)(msg + 1), memberNode->memberList.data(), listSize);



    emulNet->ENsend(&memberNode->addr, &target, (char *)msg, msgSize);
    free(msg);

    

 

}

//merge membership list recieved from other node
void MP1Node::mergeMembershipList(char *data, int size) {
    int count = size / sizeof(MemberListEntry);
    MemberListEntry *entries = (MemberListEntry *)data;

    for (int i = 0; i < count; i++) {
        Address addr = getNodeAddress(entries[i].id, entries[i].port);
        addOrUpdateMember(addr, entries[i].heartbeat, par->getcurrtime());
    }
}

//send gossip message to single target node
void MP1Node::sendGossip(Address target) {
    size_t listSize = memberNode->memberList.size() * sizeof(MemberListEntry);
    size_t msgSize = sizeof(MessageHdr) + listSize;

    MessageHdr *msg = (MessageHdr *)malloc(msgSize);
    msg->msgType = GOSSIP;
    memcpy((char *)(msg + 1), memberNode->memberList.data(), listSize);
    emulNet->ENsend(&memberNode->addr, &target, (char *)msg, msgSize);
    free(msg);
    
}




/**
* FUNCTION NAME: nodeLoopOps
*
* DESCRIPTION: Check if any node hasn't responded within a timeout period and then delete
*                 the nodes
*                 Propagate your membership list
*/
void MP1Node::nodeLoopOps() {
    
    /*
     * Your code goes here
     */

    long currTime = par->getcurrtime();

    //incrementing own heartbeat and update in member list
    memberNode->heartbeat++;
    int selfId;
    short selfPort;
    memcpy(&selfId, &memberNode->addr.addr[0], sizeof(int));
    memcpy(&selfPort, &memberNode->addr.addr[4], sizeof(short));

    for (auto &m : memberNode->memberList) {
        if(m.id == selfId && m.port == selfPort) {
            m.heartbeat = memberNode->heartbeat;
            m.timestamp = currTime;
            break;
        }
    }

    //failure detection and cleanup
    for (auto it = memberNode->memberList.begin(); it != memberNode->memberList.end();) {

        int id = it->id;
        short port = it->port;

        //skipped self
        if (id == selfId && port == selfPort) {
            ++it;
            continue;
        }

        long time_difference = currTime - it->timestamp;

        //increment timeout counter based on time since last heardd from
        //timeout_map[{id, port}] += currTime - lastHeard;

        //removing node if not heard from TREMOVE
        //if(timeout_map[{id,port}] >= TREMOVE) {
        if(time_difference > TREMOVE) {
            //Address nodeAddr;
            Address nodeAddr = getNodeAddress(id,port);
            
            //memcpy(&nodeAddr.addr[0], &id, sizeof(int));
            //memcpy(&nodeAddr.addr[4], &port, sizeof(short));

            log->logNodeRemove(&memberNode->addr, &nodeAddr);
            it = memberNode->memberList.erase(it);
            continue;
            
            //timeout_map.erase({id,port});
        }
        else {
            //timeout_map[{id,port}] = time_difference;   //added this line 
            ++it;
        }
    }

    #ifdef DEBUGLOG
            log->LOG(&memberNode->addr, "DBG after removals in nodeLoopOps:");
            debugPrintMemberList("after removals");
            #endif

    gossipMembershipList();


    #ifdef DEBUGLOG
            log->LOG(&memberNode->addr, "DBG after gossip:");
            debugPrintMemberList("after gossip");
            #endif
 
}

/**
 * FUNCTION NAME: isNullAddress
 *
 * DESCRIPTION: Function checks if the address is NULL
 */
int MP1Node::isNullAddress(Address *addr) {
	return (memcmp(addr->addr, NULLADDR, 6) == 0 ? 1 : 0);
}

/**
 * FUNCTION NAME: getJoinAddress
 *
 * DESCRIPTION: Returns the Address of the coordinator
 */
Address MP1Node::getJoinAddress() {
    Address joinaddr;

    *(int *)(&joinaddr.addr) = 1;
    *(short *)(&joinaddr.addr[4]) = 0;

    return joinaddr;
}

/**
 * FUNCTION NAME: initMemberListTable
 *
 * DESCRIPTION: Initialize the membership list
 */
void MP1Node::initMemberListTable(Member *memberNode) {
	memberNode->memberList.clear();
}

/**
 * FUNCTION NAME: printAddress
 *
 * DESCRIPTION: Print the Address
 */
void MP1Node::printAddress(Address *addr)
{
    printf("%d.%d.%d.%d:%d \n",  addr->addr[0],addr->addr[1],addr->addr[2],
                                                       addr->addr[3], *(short*)&addr->addr[4]) ;    
}



//DEBUGGG HELPER FUNC

void MP1Node::debugPrintMemberList(const char *tag) {
#ifdef DEBUGLOG
    char s[512];
    sprintf(s, "DBG %s: memberList size=%zu", tag, memberNode->memberList.size());
    log->LOG(&memberNode->addr, s);
    for (auto &m : memberNode->memberList) {
        sprintf(s, "DBG  member %d:%d hb=%ld ts=%ld",
        m.id, m.port, m.heartbeat, m.timestamp);
        log->LOG(&memberNode->addr, s);
    }

#endif
}
