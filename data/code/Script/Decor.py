#-*- coding:utf-8 -*-

from __future__ import division
from functools import wraps
import logging
logging.getLogger().setLevel(logging.INFO)


class Log:
    @classmethod
    def loginfo(self, func):
        @wraps(func)
        def log(*args, **kwargs):
            try:
                logging.info(u"当前运行方法：%s", func.__name__)
                return func(*args, **kwargs)
            except Exception as e:
                logging.warn(u"方法:%s 错误:%s", func.__name__, e)

        return log
