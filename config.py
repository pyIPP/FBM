import os

sharesRoot = '/shares/departments/AUG/users/%s' %os.getenv('USER')
udbDir = '%s/udb' %sharesRoot
tr_clientDir = '%s/tr_client/AUGD' %sharesRoot
