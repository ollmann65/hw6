//
// Wrapper around ht.h to make all private/protected functions public
//

#ifndef CS104_HW6_TEST_SUITE_PUBLICIFIED_HT_H
#define CS104_HW6_TEST_SUITE_PUBLICIFIED_HT_H

#define private public
#define protected public
#include <ht.h>
#undef private
#undef public

#endif //CS104_HW6_TEST_SUITE_PUBLICIFIED_HT_H
