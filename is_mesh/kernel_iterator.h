
#pragma once

#include <iterator>
#include <memory>
#include <vector>
#include <cassert>

#include <is_mesh/kernel.h>

namespace is_mesh
{
    
    /*
     * Forward declaration
     template<
     typename value_type
     , typename key_type
     , typename allocator_type
     > class kernel;
     */
    
    /**
     * An iterator class used by the is_mesh::kernel.
     * The kernel wraps a kernel handle and uses indirect access to the kernel.
     * Iterators are valid through out their existence - only if an element pointed
     * to by an iterator is it invalidated, but if the delete is undone the iterator
     * will be valid again after the undo operation.
     *
     * The garbage collect routine in the kernel might change the ordering of
     * elements in the kernel and provide different results is one iterates through
     * the kernel before and after a garbage collect.
     * The undo operation is guaranteed to keep the ordering before and after a mark/undo
     * pair of operations.
     *
     * @see kernel.
     */
    template <class key_t_>
    class kernel_iterator
    {
    private:
        typedef typename key_t_::kernel_element    element_type;
        //    typedef typename element_type::element_pointer  element_pointer;
    public:
        typedef          key_t_                                       kernel_type;
        typedef          kernel_iterator<kernel_type>                 iterator;
        typedef typename kernel_type::kernel_element                  kernel_element;
        typedef typename kernel_type::allocator_type::difference_type difference_type;
        typedef typename element_type::value_type                     value_type;
        typedef typename element_type::key_type                       key_type;
        typedef          value_type *                                 pointer;
        typedef          value_type &                                 reference;
        typedef          std::input_iterator_tag                      iterator_category;
        
    private:
        
        key_type         m_key;
        kernel_type*     m_kernel;
        value_type*      m_value;   //temporary storage for value when returning it
        
    public:
        /**
         * The only constructor.
         * Creates a kernel iterator. Should only be created from the kernel.
         */
        kernel_iterator(kernel_type const * const kernel, key_type const & key) : m_key(key)
        //kernel_iterator(kernel_type const * const kernel, key_type & key) : m_key(key)
        {
            //hack.. can this be avoided???
            m_kernel = (kernel_type*) &*kernel;
        }
        
    public:
        /**
         * Converts the iterator to a handle or key.
         */
        key_type     key()    const { return m_key; }
        /**
         * Returns a pointer to the kernel that the iterator is bound to.
         */
        kernel_type* kernel() const { return m_kernel; }
        
        /**
         * Compares two iterators for equality.
         */
        friend bool operator==(const iterator& i, const iterator& j)
        {
            return (i.m_key == j.m_key) && (i.m_kernel == j.m_kernel);
        }
        /**
         * Compares two iteators for inequality.
         */
        friend bool operator!=(const iterator& i, const iterator& j)
        {
            return (i.m_key != j.m_key) || (i.m_kernel != j.m_kernel);
        }
        
        /**
         * The member access operator.
         *
         * @return Pointer to the element contained within the kernel cell.
         */
        pointer operator->()
        {
            assert(m_kernel->lookup(m_key).state == element_type::VALID);
            m_value = &m_kernel->lookup(m_key).value;
            return m_value;
        }
        
        /**
         * The dereference operator.
         *
         * @return The element that is contained within the kernel cell.
         */
        reference operator*()
        {
            assert(m_kernel->lookup(m_key).state == element_type::VALID);
            m_value = &m_kernel->lookup(m_key).value;
            return *m_value;
        }
        
        /**
         * The pre-increment operator.
         *
         * @return The next iterator in the kernel sequence.
         */
        iterator& operator++()
        {
            kernel_element& cur = m_kernel->lookup(m_key);
            m_key = cur.next;
            return *this;
        }
        
        /**
         * The post-increment operator.
         *
         * @return The iterator
         */
        iterator operator++(int)
        {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }
    };
}
