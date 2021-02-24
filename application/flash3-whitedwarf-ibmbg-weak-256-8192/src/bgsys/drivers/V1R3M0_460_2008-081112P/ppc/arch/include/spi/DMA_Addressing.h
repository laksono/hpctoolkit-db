/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* (C) Copyright IBM Corp.  2008, 2008                              */
/* IBM CPL License                                                  */
/*                                                                  */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */

#ifndef	_DMA_ADDRESSING_H_ /* Prevent multiple inclusion */
#define	_DMA_ADDRESSING_H_


/*!
 * \file spi/DMA_Addressing.h
 * 
 * \brief DMA SPI Addressing Definitions and Inline Functions
 *
 * This include file contains definitions and inline functions that are used to
 * establish memory addresses and offsets used by the DMA.  
 *
 * Functions are provided to initialize addressing, obtain a handle representing
 * the local process' addressing, and resolve a virtual address associated with
 * a particular handle.
 *
 */


#include <common/namespace.h>


__BEGIN_DECLS


/*!
 * \brief __INLINE__ definition
 * 
 * Option 1:
 * Make all functions be "static inline":
 * - They are inlined if the compiler can do it
 * - If the compiler does not inline it, a single copy of the function is
 *   placed in the translation unit (eg. xxx.c)for use within that unit.
 *   The function is not externalized for use by another unit...we want this
 *   so we don't end up with multiple units exporting the same function,
 *   which would result in linker errors.
 *
 * Option 2:
 * A GNU C model: Use "extern inline" in a common header (this one) and provide
 * a definition in a .c file somewhere, perhaps using macros to ensure that the
 * same code is used in each case. For instance, in the header file:
 *
 * \verbatim
   #ifndef INLINE
   # define INLINE extern inline
   #endif
   INLINE int max(int a, int b) {
     return a > b ? a : b;
   }
   \endverbatim
 *
 * ...and in exactly one source file (in runtime/SPI), that is included in a
 * library...
 *
 * \verbatim
   #define INLINE
   #include "header.h"
   \endverbatim
 * 
 * This allows inlining, where possible, but when not possible, only one 
 * instance of the function is in storage (in the library).
 */
#ifndef __INLINE__
#define __INLINE__ extern inline
#endif


#include <spi/kernel_interface.h>   /* For Kernel_ functions          */
#include <spi/DMA_Assert.h>
#include <spi/DMA_Counter.h>
#include <spi/DMA_InjFifo.h>
#include <spi/DMA_RecFifo.h>


/*!
 * \brief DMA Addressing Handle
 * 
 * A handle, that is opaque to the DMA SPI user, representing the addressing
 * information for a process.
 *
 * \note Internally, it is the offset within the Global Application Segment Sets 
 * of the application segment set that matches this process' application segment
 * set, relative to the start of the first global application segment set.
 * It is intended that this offset be passed to other nodes along with virtual
 * addresses.  The other nodes can pass this offset into the address resolution
 * function ( DMA_Addressing_resolve() ) so it can use the proper set of 
 * application segments to perform the resolution.
 *
 * \see DMA_AddressingGetHandle
 */
typedef unsigned DMA_AddressingHandle_t;


/*!
 * \brief Pointer to Local Barrier Function Used By DMA_AddressingInitPhase2()
 *
 * This barrier function must be coded such that it performs a barrier
 * among the processes on our physical node (ie. a local process barrier).
 * When all processes on our physical node have called this function,
 * this function returns.
 */
typedef void (*DMA_AddressingInitLocalBarrierFcn_t)();


/*!
 * \brief Pointer to Global Barrier Function Used By DMA_AddressingInitPhase2()
 *
 * This barrier function must be coded such that it performs a barrier
 * among the physical nodes in the partition.
 */
typedef void (*DMA_AddressingInitGlobalBarrierFcn_t)();


#include <spi/DMA_AddressingInternals.h>


/*!
 * \brief Initialize DMA Addressing Phase 1
 *
 * Phase 1 DMA Addressing Init.  Must be called before any DMA resources are
 * allocated (such as counters and fifos).
 *
 * \retval  0            Success
 * \retval  errorNumber  Failure
 */
int DMA_AddressingInitPhase1(void);


/*!
 * \brief Initialize DMA Addressing Phase 2
 *
 * Phase 2 DMA addressing init.  Must be called after DMA resources (such as
 * counters and fifos) have been allocated and initialized on all nodes
 * in the partition.  This is a collective operation, so all nodes must 
 * participate.
 *
 * \param[in]  localBarrierFn  Function pointer that will perform a local 
 *                             process barrier, allocated using
 *                             LockBox_AllocateProcessBarrier().  This is used
 *                             to barrier among the virtual nodes that are
 *                             local to this physical node.
 * \param[in]  localBarrierParam  Opaque pointer to be passed to the
 *                                localBarrierFn.
 * \param[in]  globalBarrierFn Function pointer that will perform a global
 *                             barrier across all physical nodes in the
 *                             partition.
 * \param[in]  globalBarrierParam  Opaque pointer to be passed to the
 *                                 globalBarrierFn.
 * \param[in]  injFifoGroup  Pointer to an initialized injection fifo group.
 * \param[in]  recFifoGroup  Pointer to an initialized reception fifo group.
 * \param[in]  counterGroup  Pointer to an initialized injection counter group.
 *
 * \retval  0            Success
 * \retval  errorNumber  Failure
 *
 * \see LockBox_AllocateProcessBarrier
 */
int DMA_AddressingInitPhase2(DMA_AddressingInitLocalBarrierFcn_t   localBarrier,
			     void                                 *localBarrierParam,
			     DMA_AddressingInitGlobalBarrierFcn_t  globalBarrier,
			     void                                 *globalBarrierParam,
			     DMA_InjFifoGroup_t *injFifoGroup,
			     DMA_RecFifoGroup_t *recFifoGroup,
			     DMA_CounterGroup_t *counterGroup);


/*!
 *
 * \brief Get DMA Addressing Handle
 *
 * Return the DMA addressing handle associated with the local process.
 *
 * \returns  handle  The addressing handle
 *
 * \pre DMA_AddressingInit was successfully called.
 */
__INLINE__ DMA_AddressingHandle_t DMA_AddressingGetHandle(void)
{
  SPI_assert ( DMA_AddressingPhase2HasBeenInitialized );

  return ( DMA_AddressingLocalAppSegmentSetInfo.addressingHandle );
}


/*!
 * \brief Get Virtual Addresses for the Min and Max Physical Addresses 
 *        for User Space
 *
 * Return the virtual addresses associated with the min and max physical 
 * addresses allowed for user space for the calling process.
 *
 * \param[out]  va_min  Pointer to a pointer.  Upon return, the pointer is 
 *                      set to the virtual address associated with the
 *                      minimum physical address allowed for user space.
 * \param[out]  va_max  Pointer to a pointer.  Upon return, the pointer is 
 *                      set to the virtual address associated with the
 *                      maximum physical address allowed for user space.
 *
 * If the DMA_AddressingInit function has not been successfully called yet,
 * a value of 0 for the min and 0xFFFFFFFF max is returned.
 *
 */
__INLINE__ void DMA_AddressingGetMinMaxVa(void ** va_min,
				          void ** va_max)
{
  if ( DMA_AddressingPhase1HasBeenInitialized )
  {
    *va_min = (void*)DMA_AddressingLocalAppSegmentSetInfo.va_min;
    *va_max = (void*)DMA_AddressingLocalAppSegmentSetInfo.va_max;
  }
  else
  {
    *va_min = (void*)0;
    *va_max = (void*)0xFFFFFFFF;
  }

  return;
}


/*!
 * \brief DMA Addressing Resolve
 *
 * Given a virtual address, get the offset from the base address associated with
 * a DMA counter.
 *
 * \param[in]  c_sw     Pointer to the software counter structure.  This is only
 *                      needed when the handle is for the local process (not for
 *                      a remote process).
 * \param[in]  va       Virtual address whose offset from the counter's base is
 *                      to be returned.
 * \param[in]  length   The number of bytes in the buffer pointed to by va.
 * \param[in]  handle   Pointer to the addressing handle associated with the
 *                      process that owns the virtual address.  If the virtual
 *                      address is for the calling process, NULL may be
 *                      specified, indicating to use the addressing information
 *                      associated with the calling process.
 *
 * \retval  offset   The offset of the va from the counter's base.
 *
 * \pre DMA_AddressingInit was successfully called.
 *
 * \note It is assumed that if the handle is not for the calling process, then
 *       the remote process' counter's base address (used in calculating the 
 *       offset) is set to the smallest physical address accessible from user 
 *       space on the process associated with the handle.  This base address 
 *       is in the global app segment set array, so c_sw is not used in this
 *       case.
 *
 */
__INLINE__ unsigned int DMA_AddressingResolve(
					      const DMA_Counter_t    *c_sw,
					      void                   *va,
					      unsigned int            length,
					      DMA_AddressingHandle_t *handlePtr
					     )
{
  SPI_assert( c_sw != NULL );
  SPI_assert( va   != NULL );
  SPI_assert_debug( DMA_AddressingPhase2HasBeenInitialized );

  DMA_AddressingHandle_t                   handle;
  DMA_AddressingLocalAppSegmentSetInfo_t  *localInfo  = &DMA_AddressingLocalAppSegmentSetInfo;
  SPI_assert( localInfo );
  DMA_AddressingGlobalAppSegmentSetInfo_t *globalInfo = &DMA_AddressingGlobalAppSegmentSetInfo;
  DMA_AddressingAppSegmentSet_t           *appSegmentSet;
  DMA_AddressingAppSegment_t              *appSegmentArray;
  uint32_t                 numAppSegments;
  uint32_t                 i;
  uint32_t                 segmentVaBase=0;
  uint32_t                 offset;
  uint32_t                 counterPaBase;

  /* Determine which app segment set to use:
   * Local handle when handlePtr is NULL or the provided handle.
   */
  if ( handlePtr == NULL)
    handle = localInfo->addressingHandle;
  else
    handle = *handlePtr;

  SPI_assert(handle <= globalInfo->lastAppSegmentSetOffset);
  appSegmentSet = (DMA_AddressingAppSegmentSet_t*)
                    ((char*)globalInfo->appSegmentSets + handle);

  /* Determine which application segment the virtual address is in. */

  numAppSegments  = appSegmentSet->numAppSegments;
  appSegmentArray = (DMA_AddressingAppSegment_t*)
                       ((char*)globalInfo->appSegmentSets + 
			appSegmentSet->appSegmentOffset);

  for (i=0; i<numAppSegments; i++)
  {
    segmentVaBase = appSegmentArray[i].va_base;
    if ( ((uint32_t)va >= segmentVaBase) &&
	 ((uint32_t)va - segmentVaBase < appSegmentArray[i].length) )
      break;
  }

  SPI_assert(i < numAppSegments);

  /*
   * Make sure buffer fits in app segment.
   */
  if ( ( (uint32_t)va + length - 1 ) > appSegmentArray[i].va_max )
  {
    printf("DMA_CounterGetOffsetFromBase: Buffer 0x%08x of length %d is out of bounds.  Check length.\n",
	   (unsigned)va, length);
    SPI_assert(0);
  }
  
  /* 
   * If the addressing handle being used is our process' handle, use the offset
   * from the specified counter's base to calculate the DMA offset.
   * Otherwise, assume the counter base is the smallest physical address 
   * accessible from user space for the process whose handle is being used and
   * use that.
   */
  if ( handle == localInfo->addressingHandle )
    counterPaBase = c_sw->pa_base;
  else
    counterPaBase = appSegmentSet->minPaAccessibleFromUserMode;

  /*
   * If the base physical address of the application segment found above is
   * greater than or equal to the counter's base physical address (typical
   * case), proceed with the calculation based on that.
   *
   * Otherwise, use a slightly different calculation (see else leg).
   */
  if ( appSegmentArray[i].pa_base >= counterPaBase )
  {
    /* 
     * Calculate the offset from the counter base:
     * - offset from app segment's virtual address base (va - segmentVaBase) + 
     * - segment's physical base (shifted) - counter's base (shifted) * 16
     */
    offset = 
      ((unsigned)va - segmentVaBase) +
      ( (appSegmentArray[i].pa_base - counterPaBase) << 4 );
    
/*     printf("GetOffsetFromBase:  va=0x%08x, length=%d, offset=0x%08x, index=%d, segmentVaBase=0x%08x, appSegmentArrayPaBase=0x%08x, counterBase=0x%08x\n",(unsigned)va, length, offset, i, */
/* 	   segmentVaBase, appSegmentArray[i].pa_base, counterPaBase); */
/*     fflush(stdout); */
  }
  /*
   * Handle the case where the counter's base exceeds the app segment's base.
   * This occurs when the counter's base is set to the base of the buffer
   * rather than the min base of all the app segments.  This can only happen
   * when the handle is our process' handle (not a handle for a remote node).
   */
  else
  {
    offset = 
      ((unsigned)va - segmentVaBase) -
      ( (counterPaBase - appSegmentArray[i].pa_base) << 4 );
    
/*     printf("GetOffsetFromBase2:  va=0x%08x, length=%d, offset=0x%08x, index=%d, segmentVaBase=0x%08x, appSegmentArrayPaBase=0x%08x, counterBase=0x%08x\n",(unsigned)va, length, offset, i, */
/* 	   segmentVaBase, appSegmentArray[i].pa_base, counterPaBase); */
/*     fflush(stdout); */
  }

  return ( offset );
}


/*!
 * \brief DMA Addressing Resolve By Id
 *
 * Given a virtual address, get the offset from the base address associated with
 * a DMA counter, using the DMA counter's id.
 *
 * \param[in]  cg_ptr      Pointer to the structure previously filled in when
 *                         the counter was allocated.
 * \param[in]  counter_id  Identifier of the counter
 *                         (0 to DMA_NUM_COUNTERS_PER_GROUP-1).
 * \param[in]  va       Virtual address whose offset from the counter's base is
 *                      to be returned.
 * \param[in]  length   The number of bytes in the buffer pointed to by va.
 * \param[in]  handle   Pointer to the addressing handle associated with the
 *                      process that owns the virtual address.  If the virtual
 *                      address is for the calling process, NULL may be
 *                      specified, indicating to use the addressing information
 *                      associated with the calling process.
 *
 * \retval  offset   The offset of the va from the counter's base.
 *
 * \pre DMA_AddressingInit was successfully called.
 *
 * \note It is assumed that if the handle is not for the calling process, then
 *       the counter's base address (used in calculating the offset) is set to
 *       the smallest physical address accessible from user space on the process
 *       associated with the handle.
 *
 */
__INLINE__ unsigned int DMA_AddressingResolveById(
				const DMA_CounterGroup_t *cg_ptr,
				int                       counter_id,
				void                     *va,
				unsigned int              length,
				DMA_AddressingHandle_t   *handle
					     )
{
  SPI_assert( (counter_id >= 0) && (counter_id < DMA_NUM_COUNTERS_PER_GROUP) );
  SPI_assert( cg_ptr != NULL );
  SPI_assert( (cg_ptr->permissions[DMA_COUNTER_GROUP_WORD_ID(counter_id)] & 
	   _BN(DMA_COUNTER_GROUP_WORD_BIT_ID(counter_id))) != 0 );

/*   printf("Getting offset from counter %d for core %d\n",counter_id,coreNum); */
  return( DMA_AddressingResolve ( &cg_ptr->counter[counter_id],
				  va,
				  length,
				  handle ) );
}


__END_DECLS


#endif 
