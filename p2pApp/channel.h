#ifndef CHANNEL_H
#define CHANNEL_H

#include <pv/pvAccess.h>

#include "chancache.h"

struct GWChannel : public epics::pvAccess::Channel,
        std::tr1::enable_shared_from_this<GWChannel>
{
    ChannelCacheEntry::shared_pointer entry;
    epics::pvAccess::ChannelRequester::shared_pointer requester;

    GWChannel(ChannelCacheEntry::shared_pointer e,
              epics::pvAccess::ChannelRequester::shared_pointer);
    virtual ~GWChannel();


    // for Requester
    virtual std::string getRequesterName();
    virtual void message(std::string const & message, epics::pvData::MessageType messageType);

    // for Destroyable
    virtual void destroy(){}

    // for Channel
    virtual std::tr1::shared_ptr<epics::pvAccess::ChannelProvider> getProvider();
    virtual std::string getRemoteAddress();

    virtual ConnectionState getConnectionState();
    virtual std::string getChannelName();
    virtual std::tr1::shared_ptr<epics::pvAccess::ChannelRequester> getChannelRequester();
    virtual bool isConnected();

    virtual void getField(epics::pvAccess::GetFieldRequester::shared_pointer const & requester,
                          std::string const & subField);
    virtual epics::pvAccess::AccessRights getAccessRights(epics::pvData::PVField::shared_pointer const & pvField);
    virtual epics::pvAccess::ChannelProcess::shared_pointer createChannelProcess(
            epics::pvAccess::ChannelProcessRequester::shared_pointer const & channelProcessRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvAccess::ChannelGet::shared_pointer createChannelGet(
            epics::pvAccess::ChannelGetRequester::shared_pointer const & channelGetRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvAccess::ChannelPut::shared_pointer createChannelPut(
            epics::pvAccess::ChannelPutRequester::shared_pointer const & channelPutRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvAccess::ChannelPutGet::shared_pointer createChannelPutGet(
            epics::pvAccess::ChannelPutGetRequester::shared_pointer const & channelPutGetRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvAccess::ChannelRPC::shared_pointer createChannelRPC(
            epics::pvAccess::ChannelRPCRequester::shared_pointer const & channelRPCRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvData::Monitor::shared_pointer createMonitor(
            epics::pvData::MonitorRequester::shared_pointer const & monitorRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);
    virtual epics::pvAccess::ChannelArray::shared_pointer createChannelArray(
            epics::pvAccess::ChannelArrayRequester::shared_pointer const & channelArrayRequester,
            epics::pvData::PVStructure::shared_pointer const & pvRequest);

    virtual void printInfo();
    virtual void printInfo(std::ostream& out);

};

#endif // CHANNEL_H
